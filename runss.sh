#!/bin/bash

#for i in {0..9}
#do 
#Umin=$(bc<<<"1.0")	
#Umid=$(bc<<<"3.0")	
#Uhigh=$(bc<<<"5.0")	
#Utop=$(bc<<<"7.0")	
for k in {0..9}
do
	for j in {0..5}
	do
			U=$(bc<<<"$j+2")	
			ms=$(bc<<<"30")
			T=$(bc<<<"$k")
			o=$(bc<<<"4000+($k)*100+$j")
			echo $o
#			./a.out -lvl 6 -ms $ms -U1 $U -U2 1 -U3 $U -U4 1 -U5 $U -U6 1 -T1 $T -T2 $T -T3 $T -T4 $T -T5 $T -T6 $T -gamma 1 -o $o -fullT 500 -dT 0.00001 -pot 2 -f $f -wf 1000
			for l in {21..302..10}
			do
				python check_markov.py -o $o -l $l -rc 1
			done
	done
done
#done






