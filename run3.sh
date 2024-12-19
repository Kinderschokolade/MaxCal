#!/bin/bash

for j in {5..5}
do
        for i in {0..9}
        do
		for t in {0..0}
		do
                #echo $o
#               T=$(bc<<<"1 + $j")
                f=$(bc<<<"0")
                U1=$(bc<<<"0.5")
                U2=$(bc<<<"6.5")
                U3=$(bc<<<"5")
                U4=$(bc<<<"5.3")
                U5=$(bc<<<"4")
                U6=$(bc<<<"4.7")
                U7=$(bc<<<"1")
                U8=$(bc<<<"1.1")
                T=$(bc<<<"$i * 0.2 + 0.2")
                dt=$(bc<<<"0.00001")
                gamma=$(bc<<<"1")

                fullT=$(bc<<<"50000")

                o=$(bc<<<"8000+ $i*100+$j*10")
		echo $o 
	     	./main_floc -lag 100 200 400 1000  -lvl 8 -U1 $U1  -U2 $U2 -U3 $U3 -U4 $U4 -U5 $U5 -U6 $U6 -U7 $U7 -U8 $U8 -T $T -gamma $gamma -seed $o -pot 1 -f $f -wf 0 -fullT $fullT -o $o -ms 30 -dT $dt     

#		python check_markov.py -o $o -l 100 -rc 1
		#python check_markov.py -o $o -l 200 -rc 1
		done
	done
done


