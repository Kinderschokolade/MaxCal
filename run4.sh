#!/bin/bash


for i in {0..9}
do
        for j in {0..0}
        do
		for t in {0..0}
		do
                #echo $o
#               T=$(bc<<<"1 + $j")
                f=$(bc<<<"$i*2")
                U1=$(bc<<<"12")
                U2=$(bc<<<"4")
                U3=$(bc<<<"8")
                U4=$(bc<<<"2")
                U5=$(bc<<<"8")
                U6=$(bc<<<"4")
                T=$(bc<<<"2")
                dt=$(bc<<<"0.00001")
                gamma=$(bc<<<"1")

                fullT=$(bc<<<"40000")

                o=$(bc<<<"2012+ $i*100+$j*10")
		echo $o 
	    	./main -lag 100 200 400 600 800 1000 1400 1800  -lvl 6 -U1 $U1 -U2 $U2 -U3 $U3 -U4 $U4 -U5 $U5 -U6 $U6 -T $T -gamma $gamma -pot 1 -f $f -wf 0 -fullT $fullT -o $o -ms 60 -dT $dt     
#		python check_markov.py -o $o -l 100 -rc 1
		python check_markov.py -o $o -l 200 -rc 1
		python check_markov.py -o $o -l 400 -rc 1
		done
	done
done


for i in {0..9}
do
        for j in {0..0}
        do
		for t in {0..0}
		do
                #echo $o
#               T=$(bc<<<"1 + $j")
                f=$(bc<<<"$i*2")
                U1=$(bc<<<"12")
                U2=$(bc<<<"4")
                U3=$(bc<<<"8")
                U4=$(bc<<<"2")
                U5=$(bc<<<"8")
                U6=$(bc<<<"4")
                T=$(bc<<<"4")
                dt=$(bc<<<"0.00001")
                gamma=$(bc<<<"1")

                fullT=$(bc<<<"40000")

                o=$(bc<<<"2013+ $i*100+$j*10")
		echo $o 
	    	./main -lag 100 200 400 600 800 1000 1400 1800  -lvl 6 -U1 $U1 -U2 $U2 -U3 $U3 -U4 $U4 -U5 $U5 -U6 $U6 -T $T -gamma $gamma -pot 1 -f $f -wf 0 -fullT $fullT -o $o -ms 60 -dT $dt     
#		python check_markov.py -o $o -l 100 -rc 1
		python check_markov.py -o $o -l 200 -rc 1
		python check_markov.py -o $o -l 400 -rc 1
		done
	done
done

