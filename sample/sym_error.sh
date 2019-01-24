#!/bin/bash

# test scaling with steps and error (energy, time)
source ~/.uftools.sh

#p2.test
tt=p2
sbase=0.49672941329
#p3.test
tt=p3
sbase=0.000015581792845073965 # use formula in Symbolic_calculator.ipynb


exec="./hierarchy $tt.test -r 2"
#kepc='../../P3TARC/test/keplersolvertest 2'

nn=4
s=(`powseq 0.5 -n $nn`)
step=(`powseq -a $sbase 0.5 -n $nn`)
n=(`powseq 2 -n $nn`)

for i in `seq 0 $nn`
do
    nstep=`echo ${n[i]}'*6400'|bc`
    echo 's='${s[i]}' step='${step[i]}' nstep='$nstep
    $exec -k -6 -s ${step[i]} -n $nstep 1>sym_${tt}_w6.s${s[i]} 2>sym_${tt}_w6_err.s${s[i]}
done

#cut -d' ' -f1 sym_p2_w6.s${s[5]}>tlist
#$kepc p2.line tlist 1>kepler_p2 2>kepler_p2_err

