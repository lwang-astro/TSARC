#!/bin/bash

# symplectic integrator benchmark analsys

# test scaling with steps
source ~/.uftools.sh

exec='./chain -N 5 input '
#exec='./hierarchy p3.test -r 2'

step=(`powseq 0.5 -s 10 -n 15`)
n=(`powseq 2 -n 5`)
order=(`seq 2 2 10` -6 -8)

## test scaling with order
for j in `seq 0 6`
do
    rm -f sym_step_k${order[j]}
    for i in `seq 0 4`
    do
 	echo 'Step='${step[i]}' n='${n[i]}' order='${order[j]}
 	$exec -k ${order[j]} -s ${step[i]} -n ${n[i]} 1>data.log 2>data.err
 	echo ${step[i]}' '`tail -1 data.log` >>sym_step_k${order[j]}
    done
done

# test long-term error

rm -f sym_prof
for j in `seq 0 4`
do
    for i in `seq 0 6`
    do
	echo 'order='${order[i]}' step='${step[j]}
	$exec -k ${order[i]} -s ${step[j]} -n ${n[j]}000 1>sym_data.k${order[i]}.s${step[j]} 2>sym_err.k${order[i]}.s${step[j]}
	egrep Time_profile sym_err.k${order[i]}.s${step[j]} |tail -1 |awk -v o=${order[i]} -v s=${step[j]} '{print o,s,$3,$13}' >>sym_prof
    done
done

echo 'Extrapolation integrator'
for i in `seq 2 4`
do
    echo 'sequence '$i
    $exec -s ${step[0]} -n ${n[0]}000 -m 'rational' -k $i 1>extra_data_s$i  2>extra_err_s$i
    egrep -i profile extra_err_s$i |tail -1 |awk '{print 0,0,$3,$11}' >>sym_prof
done
