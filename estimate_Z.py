import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
args = parser.parse_args();
ident = args.o


data = np.loadtxt("/data/pckr146/bause/single_particle/trajectory_"+ str(ident) + ".dat", usecols = (1,), delimiter="\t")

print(data[0], data[1])
