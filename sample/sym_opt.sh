#!/bin/bash


# test scaling with steps
source ~/.uftools.sh

step=(`powseq 0.5 -s 10 -n 15`)
n=(`powseq 2 -n 5`)

for i in `seq 0 4`
do
    echo 'step='${step[i]}
    for j in `seq 0 4`
    do
	echo 'e'$j
	./hierarchy p3.e$j -k -6 -s ${step[i]} -n ${n[i]}000 1>sym_data.s${step[i]}.e$j 2>sym_err.s${step[i]}.e$j
    done
done
