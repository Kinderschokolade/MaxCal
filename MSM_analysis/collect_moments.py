import argparse

import numpy as np
from define_flow import define_flow_direction, read_cut
from fct_markov import analyse_MSM, calc_av, calc_Sprod

parser = argparse.ArgumentParser()
parser.add_argument("-o", help="start from here")
parser.add_argument("-rc", help="rancut")
parser.add_argument("-bstart", help="xbxx")
parser.add_argument("-bstop", help="xbxx")
parser.add_argument("-cstart", help="xxcx")
parser.add_argument("-cstop", help="xxcx")
parser.add_argument("-dstart", help="xxxd")
parser.add_argument("-dstop", help="xxxd")
parser.add_argument("-out", help="out ident")
parser.add_argument("-l", help="lagtime")
args = parser.parse_args()
start = str(args.o)
bstart = int(args.bstart)
bstop = int(args.bstop)
cstart = int(args.cstart)
dstart = int(args.dstart)
cstop = int(args.cstop)
dstop = int(args.dstop)
out = str(args.out)
l = args.l
rancut = int(args.rc)


minpos = [1, 9, 16, 23]

path0 = "/data/isilon/bause/"
for line in open(path0 + "single_particle/param_" + str(start) + ".dat", "r"):
    cont = line.split()
    if cont[0] == "microstates":
        ms = int(cont[1])

temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/MFPT/moment1_" + str(start) + ".dat")
lvl = np.shape(temp)[0]
print("level", lvl)


b = bstop - bstart + 1
c = cstop - cstart + 1
d = dstop - dstart + 1


mom1 = np.zeros((lvl * lvl + 1, c * b * d))
mom2 = np.zeros((lvl * lvl + 1, c * b * d))
mom3 = np.zeros((lvl * lvl + 1, c * b * d))
pstat = np.zeros((lvl + 1, c * b * d))
Sprod = np.zeros((2, c * b * d))

ms, cut = read_cut(path0, start)
F = define_flow_direction(ms, cut)

for i in range(b):
    for j in range(c):
        for k in range(d):
            ident = str(int(start) + i * 100 + j * 10 + k)
            print(ident)

            for line in open(path0 + "single_particle/param_" + str(ident) + ".dat", "r"):
                cont = line.split()
                if cont[0] == "extf":
                    force = int(cont[1])

            T = np.loadtxt(path0 + "single_particle/" + str(l) + "/T/T_" + ident + ".dat")
            ms = np.shape(T)[0]
            for m in range(ms):
                for n in range(ms):
                    if np.isnan(T[m, n]):
                        T[m, n] = 0.0
            Ev, Evec = analyse_MSM(T)
            p = Evec[0] / sum(Evec[0])
            J = calc_av(F, T, p)
            temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/MFPT/moment1_" + ident + ".dat")
            mom1[1:, i + j * b + b * c * k] = temp.reshape(lvl * lvl)
            mom1[0, i + j * b] = ident

            temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/MFPT/moment2_" + ident + ".dat")
            mom2[1:, i + j * b + b * c * k] = temp.reshape(lvl * lvl)
            mom2[0, i + j * b] = ident

            temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/MFPT/moment3_" + ident + ".dat")
            mom3[1:, i + j * b + b * c * k] = temp.reshape(lvl * lvl)
            mom3[0, i + j * b] = ident

            Sprod[0, i + j * b] = force
            Sprod[1, i + j * b] = calc_Sprod(T)

            temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/EV/ps_" + ident + ".dat")
            pstat[0, i + j * b] = force
            for k in range(lvl):
                pstat[k + 1, i + j * b] = sum(temp[minpos[k] - rancut : minpos[k] + rancut + 1])


np.savetxt(path0 + "single_particle/" + str(l) + "/MFPT/mom3_" + out + ".dat", np.transpose(mom3))
np.savetxt(path0 + "single_particle/" + str(l) + "/MFPT/mom2_" + out + ".dat", np.transpose(mom2))
np.savetxt(path0 + "single_particle/" + str(l) + "/MFPT/mom1_" + out + ".dat", np.transpose(mom1))
np.savetxt(path0 + "single_particle/" + str(l) + "/MFPT/Sprod_" + out + ".dat", np.transpose(Sprod))

np.savetxt(path0 + "single_particle/" + str(l) + "/T/P_" + out + ".dat", np.transpose(pstat))


momc1 = np.zeros((11, c * b))
momc2 = np.zeros((11, c * b))
momc3 = np.zeros((11, c * b))

for i in range(b):
    for j in range(c):
        ident = str(int(start) + i * 100 + j * 10)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment1_" + ident + "_" + out + ".dat")
        momc1[0, i + j * b] = ident
        momc1[1, i + j * b] = int(out)
        momc1[2:, i + j * b] = temp.reshape(lvl * lvl)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment2_" + ident + "_" + out + ".dat")
        momc2[0, i + j * b] = ident
        momc2[1, i + j * b] = int(out)
        momc2[2:, i + j * b] = temp.reshape(lvl * lvl)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment3_" + ident + "_" + out + ".dat")
        momc3[0, i + j * b] = ident
        momc3[1, i + j * b] = int(out)
        momc3[2:, i + j * b] = temp.reshape(lvl * lvl)


np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom1_x-" + out + ".dat", np.transpose(momc1))
np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom2_x-" + out + ".dat", np.transpose(momc2))
np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom3_x-" + out + ".dat", np.transpose(momc3))

momc1 = np.zeros((11, c * b))
momc2 = np.zeros((11, c * b))
momc3 = np.zeros((11, c * b))

for i in range(b):
    for j in range(c):
        ident = str(int(start) + i * 100 + j * 10)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment1_" + out + "_" + ident + ".dat")
        momc1[0, i + j * b] = ident
        momc1[1, i + j * b] = int(out)
        momc1[2:, i + j * b] = temp.reshape(lvl * lvl)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment2_" + out + "_" + ident + ".dat")
        momc2[0, i + j * b] = ident
        momc2[1, i + j * b] = int(out)
        momc2[2:, i + j * b] = temp.reshape(lvl * lvl)

        temp = np.loadtxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/moment3_" + out + "_" + ident + ".dat")
        momc3[0, i + j * b] = ident
        momc3[1, i + j * b] = int(out)
        momc3[2:, i + j * b] = temp.reshape(lvl * lvl)


np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom1_" + out + "-x.dat", np.transpose(momc1))
np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom2_" + out + "-x.dat", np.transpose(momc2))
np.savetxt(path0 + "single_particle/" + str(l) + "/reweight/MFPT/mom3_" + out + "-x.dat", np.transpose(momc3))
