#!/bin/bash

# test scaling with steps
source ~/.uftools.sh

#p2.test
tt=p2
sbase=2.8099258924162904
#p3.test
tt=p3
sbase=0.140496294621

exec="./hierarchy $tt.test -r 2"
#kepc='../../P3TARC/test/keplersolvertest 2'

nn=10
s=(`powseq 0.5 -n $nn`)
step=(`powseq -a $sbase 0.5 -n $nn`)
n=(`powseq 2 -n $nn`)

for i in `seq 0 $nn`
do
    echo 's='${s[i]}' step='${step[i]}' nstep='${n[i]}0
    $exec -k -6 -s ${step[i]} -n ${n[i]}0 1>sym_${tt}_w6.s${s[i]} 2>sym_${tt}_w6_err.s${s[i]}
done

#cut -d' ' -f1 sym_p2_w6.s${s[5]}>tlist
#$kepc p2.line tlist 1>kepler_p2 2>kepler_p2_err

