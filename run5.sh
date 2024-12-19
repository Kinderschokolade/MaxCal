#!/bin/bash


for i in {0..0}
do  
	for j in {0..0}
	do
		for k in {0..9}
		do
			for p in {0..9}
			do
			for l in {80..81..10}
				do
					in=$(bc<<<"4000+$k*100+$j+$i*10")
					o=$(bc<<<"4000+$p*100+$j+$i*10")
					python reweight_MSM_J.py -o $o -i $in -rc 1 -l $l
				done
			done
			#python Fvsgamma.py -o $o -i $in -rc 1 -l 40 
		done
	done
	echo "set $i finished"
done


#for i in {0..10}
#do  
#	T=$(bc<<<"10+$i")
#	o=$(bc<<<"50+2*$i")
#	#$echo $T
#	echo $o
#	python alpha.py -T1 5 -T2 10 -T3 $T -U1 5 -U2 1 -U3 1 -gamma 1 -dT 0.0002 -fullT 1000000 -k 100 -o $o
#done
