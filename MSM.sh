#!/bin/bash

for i in {0..10}
do
	for j in {0..9}
	do
		pos=$(bc<<<"2000 + $i*200 +$j")
		echo $pos
		python minimise_fct.py -o $pos
	done
done
