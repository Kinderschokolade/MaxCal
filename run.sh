#!/bin/bash

for j in {7..7}
do
        for i in {0..0}
        do
		for t in {0..9}
		do
                #echo $o
#               T=$(bc<<<"1 + $j")
                f=$(bc<<<"$i")
                U1=$(bc<<<"6")
                U2=$(bc<<<"2")
                U3=$(bc<<<"4")
                U4=$(bc<<<"1")
                U5=$(bc<<<"4")
                U6=$(bc<<<"2")
                T=$(bc<<<"1")
                dt=$(bc<<<"0.00001")
                gamma=$(bc<<<"1")

                fullT=$(bc<<<"($t+1)*1000")

                o=$(bc<<<"1000+ $i*100+$j*10+$t")
		echo $o 
	     	./main_fast -lag 100 200 -lvl 6 -U1 $U1 -U2 $U2 -U3 $U3 -U4 $U4 -U5 $U5 -U6 $U6 -T $T -gamma $gamma -seed $o -pot 1 -f $f -wf 0 -fullT $fullT -o $o -ms 60 -dT $dt     

		#python check_markov.py -o $o -l 200 -rc 2 -db 1
		#python create_hist.py -o $o -l 100
		#python reweight_MSM_Sloc.py -i $o -o $o -l 200 -rc 1
		done
	done
done


for i in {1080..1089}; do python -W ignore  create_MSM.py -o $i -rc 1 -l 100 ; done
