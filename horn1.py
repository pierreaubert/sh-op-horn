#!/usr/bin/env python3
#                                                  -*- coding: utf-8 -*-
# A library to optimize horns
#
# Copyright (C) 2021 Pierre Aubert pierreaubert(at)yahoo(dot)fr
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys, os, os.path

import argparse
import cmath
from math import pi, sin, cos, asin, log10, log2
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt
from sympy import solve, symbols
import numpy as np
import pandas as pd
import altair as alt
from altair_saver import save as alt_save
import gmsh

import firedrake as fd
from firedrake.petsc import PETSc

sys.path.append(os.path.expanduser("/home/pierre/src/spinorama/src"))
from spinorama.graph import (
    graph_contour,
    contour_params_default,
    graph_isoband,
    isoband_params_default,
    graph_spinorama,
    graph_params_default,
)
from spinorama.load_misc import sort_angles, graph_melt
from spinorama.compute_cea2034 import compute_cea2034
from spinorama.compute_normalize import (
    normalize_mean,
    normalize_cea2034,
    normalize_graph,
)


def helmotz_complex(mesh, freq: int, conf: Dict):
    c_air = conf["c_air"]
    rho_air = conf["rho_air"]
    g_inlet = conf["g_inlet"]
    g_outlet = conf["g_outlet"]
    V = fd.FunctionSpace(mesh, "CG", 1)
    u = fd.TrialFunction(V)
    v = fd.TestFunction(V)
    kappa = fd.Constant(2 * pi * freq / c_air)
    a_prod = fd.inner(fd.grad(u), fd.grad(v)) - kappa ** 2 * fd.inner(u, v)
    # right hand side
    a_in = rho_air * 1j * kappa * fd.inner(u, v)
    a_out = (rho_air * 1j * kappa + 1 / 2) * fd.inner(u, v)
    a = a_prod * fd.dx + a_in * fd.ds(g_inlet) + a_out * fd.ds(g_outlet)
    L = -2 * rho_air * 1j * kappa * fd.conj(v) * fd.ds(g_inlet)
    # solve
    # assemble(a)
    # assemble(L)
    sol = fd.Function(V)
    fd.solve(
        a == L,
        sol,
        solver_parameters={
            "ksp_type": "gmres",
            # "mat_type": "matfree",
            # "ksp_type": "cg",
            # "ksp_monitor": None,
            # "pc_type": "python",
            # "pc_python_type": "firedrake.AssembledPC",
            # "assembled_pc_type": "ilu",
        },
    )
    return sol


def compute_H_measurements(freqs: List[int], results):
    def angle(theta: float):
        if theta == 0:
            return "On Axis"
        return "{}°".format(int(theta))

    r = 0.99
    dfs = []
    for theta in np.linspace(0, 180, 19):  # use s/19/37/ for 5 degrees increment
        dbs = []
        theta_rad = theta * pi / 180
        p_x = r * cos(theta_rad)
        p_y = r * sin(theta_rad)
        for fr, sol in zip(freqs, results):
            p_p = sol.at(p_x, p_y, dont_raise=True)
            if p_p is not None and abs(p_p) >= 0:
                if abs(p_p) == 0:
                    p_db = -120
                else:
                    p_db = 105 + 20 * log10(abs(p_p))
                # print('{:3.0f} {:+0.2f} {:+0.2f} {}'.format(theta, p_x, p_y, p_db))
                dbs.append(p_db)
            else:
                print("{} {} {} {} ERROR".format(theta, p_x, p_y, p_p))
        if theta == 0:
            dfs.append(pd.DataFrame({"Freq": freqs, angle(theta): dbs}))
        else:
            angle_p = angle(theta)
            angle_m = "-{}".format(angle_p)
            dfs.append(pd.DataFrame({angle_p: dbs}))
            if theta != 180:
                dfs.append(pd.DataFrame({angle_m: dbs}))

    df = sort_angles(pd.concat(dfs, axis=1))
    return df, graph_melt(df)


def save_H_measurements(dfs, freq_min, freq_max):
    return (
        alt.Chart(dfs)
        .mark_line()
        .transform_filter(
            alt.FieldOneOfPredicate(
                field="Measurements",
                oneOf=[
                    "On Axis",
                    "10°",
                    "20°",
                    "30°",
                    "40°",
                    "50°",
                    "60°",
                    "90°",
                    "120°",
                    "150°",
                ],
            ),
        )
        .encode(
            x=alt.X(
                "Freq:Q",
                scale=alt.Scale(
                    type="log", base=10, nice=False, domain=[freq_min, freq_max]
                ),
            ),
            y=alt.Y("dB:Q", scale=alt.Scale(type="linear", domain=[60, 110])),
            color=alt.Color("Measurements"),
        )
        .properties(width=600)
    )


def save_H_contour(df, freq_min, freq_max):
    isoband_params = isoband_params_default
    isoband_params["xmin"] = freq_min
    isoband_params["xmax"] = freq_max
    return graph_isoband(df, isoband_params)


def save_J(freqs, results, conf, freq_min, freq_max):
    g_inlet = conf["g_inlet"]
    length_in = abs(fd.assemble(results[0] / results[0] * fd.ds(g_inlet)))
    reflections = [
        abs(abs(fd.assemble(abs(s_hz) * fd.ds(g_inlet)) / length_in) - 1)
        for s_hz in results
    ]

    return (
        alt.Chart(pd.DataFrame({"Freq": freqs, "Reflection": reflections}))
        .mark_line()
        .encode(
            x=alt.X(
                "Freq",
                title="Freq (Hz)",
                scale=alt.Scale(
                    type="linear", base=10, nice=False, domain=[freq_min, freq_max]
                ),
            ),
            y=alt.Y(
                "Reflection",
                scale=alt.Scale(type="log", base=10, nice=False, domain=[0.001, 2]),
            ),
        )
    )


def save_spin(df):
    spin = compute_cea2034(df, df)
    mean = normalize_mean(spin)
    n_spin = normalize_cea2034(graph_melt(spin), mean)
    return graph_spinorama(n_spin, graph_params_default)


def save_results(freqs, results, radical):
    for f, r in zip(freqs, results):
        output = fd.File("output/{}-{:04d}.pvd".format(radical, int(f)))
        output.write(r, t=f)


def save_graphs(freqs, results, radical, conf):
    save_results(freqs, results, radical)
    df, dfs = compute_H_measurements(freqs, results)
    freq_min = conf["freq_min"]
    freq_max = conf["freq_max"]
    graph_J = save_J(freqs, results, conf, freq_min, freq_max)
    graph_H_measurements = save_H_measurements(dfs, freq_min, freq_max)
    graph_H_contour = save_H_contour(df, freq_min, freq_max)
    graph_spin = save_spin(df)
    for name, graph in [
        ("J", graph_J),
        ("H_spl", graph_H_measurements),
        ("H_contour", graph_H_contour),
        ("spin", graph_spin),
    ]:
        alt_save(graph, "output/{}_{}.png".format(radical, name))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", help="mesh file")
    parser.add_argument(
        "--freq_min", default=50, type=int, help="compute over frequencies"
    )
    parser.add_argument(
        "--freq_max", default=1050, type=int, help="compute over frequencies"
    )
    parser.add_argument("--freq_steps", default=101, type=int, help="number of steps")
    args = parser.parse_args()

    solver_config = {}
    solver_config["freq_min"] = args.freq_min
    solver_config["freq_max"] = args.freq_max
    solver_config["freq_steps"] = args.freq_steps
    # speed of sound in air (m/s)
    solver_config["c_air"] = 343.344
    #
    solver_config["rho_air"] = 1.20458
    # incoming wave amplitude
    solver_config["p_A_in"] = 1
    # outgoing wave amplitude
    solver_config["p_B_out"] = 1
    #
    solver_config["g_inlet"] = 2
    solver_config["g_outlet"] = 3

    if args.mesh is None:
        print("--mesh is required")
        sys.exit(-1)

    results = []
    freqs = np.logspace(
        log10(solver_config["freq_min"]),
        log10(solver_config["freq_max"]),
        solver_config["freq_steps"],
    )

    comm = fd.COMM_WORLD

    PETSc.Sys.Print("setting up mesh across %d processes" % fd.COMM_WORLD.size)

    mesh = fd.Mesh("meshes/{}.msh".format(args.mesh), comm=comm)

    PETSc.Sys.Print(
        "  rank {} owns {} elements and can access {} vertices".format(
            mesh.comm.rank, mesh.num_cells(), mesh.num_vertices(), comm=fd.COMM_SELF
        )
    )

    for i, f_hz in enumerate(freqs):

        s_hz = helmotz_complex(mesh, f_hz, solver_config)
        results.append(s_hz)

        if i > 0 and (i * 100 // args.freq_steps) % 10 == 0:
            PETSc.Sys.Print(
                "{}%".format(i * 100 // args.freq_steps), end="", flush=True
            )
        if i % 4 == 0:
            PETSc.Sys.Print(".", end="", flush=True)
    PETSc.Sys.Print("done")

    if comm.rank == 0:
        save_graphs(freqs, results, args.mesh, solver_config)
    sys.exit(0)
