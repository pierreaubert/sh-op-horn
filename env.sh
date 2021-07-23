#!/bin/sh

touch env.log

mkdir -p output meshes

# to please firedrack
export OMP_NUM_THREADS=1

# to please pyright
export PYTHONPATH=$HOME/src/spinorama

## SSH AGENT
## ----------------------------------------------------------------------
ssh-agent -k 2>&1 >> env.log
eval `ssh-agent`
echo $SSH_AGENT_SOCK
if ! test -f ~/.ssh/id_rsa_github; then
    echo "ERROR github key don\'t exists!"
fi

## Github keys
## ----------------------------------------------------------------------
github=$(ssh-add -l | grep github | cut -d ' ' -f 3)
if test -z $github; then
    ssh-add ~/.ssh/id_rsa_github 2>&1 >> env.log
    github=$(ssh-add -l 2>&1 | grep github | cut -d ' ' -f 3)
fi

## python virtualenv
## ----------------------------------------------------------------------
FIREDRAKE=$HOME/src/firedrake
export PYTHONPATH=$FIREDRAKE/src:$PYTHONPATH
source $FIREDRAKE/bin/activate

## summary
## ----------------------------------------------------------------------
echo 'HSO          ' $SPIN 
echo '  '$(python3 --version) $(which python3)
echo '  '$(pip3 -V)
echo '  jupyter-lab ' $(jupyter-lab --version) $(which jupyter-lab)
echo '  PYTHONPATH  ' $PYTHONPATH
echo '  github key  ' $github
