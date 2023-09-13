#!/bin/bash

source /home/jaya/softwares/.gmx_singularity.bash

mdp=$1
tpr=$2
gro=$3
top=$4
tag=$5

gmx_s grompp -f $mdp -o $tpr -c $gro -r $gro -p ${top}

gmx_d mdrun -s $tpr -rerun $gro -deffnm rerun_$tag

gmx_d energy -f rerun_${tag}.edr -o energy_${tag}.xvg
