import argparse
import getpass
import math

import numpy as np
from analyse_trajectory import (
    mfpt_trajectory_area,
    mfpt_trajectory_ms,
    mfpt_trajectory_ms_cross,
)
from define_cross import define_cross
from define_flow import define_flow_direction, read_cut, read_minpos
from msmtools.analysis import mfpt
from pyemma.msm import MSM

from fpt_ana import first_passage

# uses -o count matrices (for different lagtime) from isilon/single_particle/count_????.dat to find
# time scale. Use to identify markovian region
# choose -l to fix lagtime and perform principal component analysis.
# Imag parts are calculated too


def flow_direction(pl, mi, ms):
    r = np.zeros((ms, ms))
    for i in range(ms):
        for j in range(ms):
            if pl[i, j] > mi[i, j]:
                r[i, j] = +1
            else:
                r[i, j] = -1
    for i in range(ms):
        r[i, i] = 0

    return r


parser = argparse.ArgumentParser()
parser.add_argument("-o", help="ident")
parser.add_argument("-l", help="write T of specific line")
parser.add_argument("-n", help="number of (slowest) eigenvectors", default=3)
parser.add_argument("-db", help="set 1 for detbal check", default=0)
parser.add_argument("-rc", help="rancut", default=1)
args = parser.parse_args()
ident = args.o
Enum = int(args.n)
wline = args.l
db_flag = int(args.db)
rancut = int(args.rc)

user = getpass.getuser()

if user == "mbause":
    path0 = "/u/mbause/data/"
    path1 = "/u/mbause/data/"
else:
    path0 = "/data/isilon/bause/"
    path1 = "/data/pckr194/bause/"


if wline != None:
    wline = int(wline)

string = " "
ms = 0
dT = 0
for line in open(path0 + "single_particle/param_" + str(ident) + ".dat", "r"):
    cont = line.split()
    if cont[0] == "dT":
        dT = float(cont[1])
    if cont[0] == "microstates":
        ms = int(cont[1])


c_s = np.zeros(ms * ms)
p_s = np.zeros(ms)

F = open(path0 + "single_particle/count_" + str(ident) + ".dat", "r")
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime = [float(i) for i in lagtime]

cols = len(lagtime)

data = np.loadtxt(path0 + "single_particle/count_" + str(ident) + ".dat", dtype="i", usecols=(range(1, cols + 1)))


ms, minpos = read_minpos(path0, ident)
ms, cut = read_cut(path0, ident)

cpl = data[0 : ms * ms]
cmi = data[ms * ms : 2 * ms * ms]
p = data[2 * ms * ms : 2 * ms * ms + ms]

L = np.zeros((ms, ms))
L = np.zeros((ms, ms))

lt = np.zeros((ms + 1, cols))

EV_arr = np.zeros((ms + 1, cols))
EV_comp = np.zeros((ms, 2))
Evec_comp = np.zeros((ms, 2))
Evecr_comp = np.zeros((ms, 2))
Evec_arr = np.zeros((cols, ms, Enum))
Evecr_arr = np.zeros((cols, ms, Enum))


# print(c.shape[1])

T_check = np.zeros((ms, ms))


c = cpl + cmi  # add all data for now..

# mfpt_ar = np.zeros((len(minpos),len(minpos)))
# mfpt_meas = np.zeros((len(minpos),len(minpos)))
# mfpt_err = np.zeros((len(minpos),len(minpos)))
# mfpt_steps = np.zeros((len(minpos),len(minpos)))

lvl = len(minpos)

mfpt_ar = np.zeros((lvl, lvl))
mfpt_arx = np.zeros((ms, ms))
mfpt_meas = np.zeros((ms, ms))
mfpt_err = np.zeros((ms, ms))
mfpt_steps = np.zeros((ms, ms))


blocks = 10
tra_len = 1000000

# r = flow_direction(cpl[:,wline].reshape((ms,ms)), cmi[:,wline].reshape((ms,ms)), ms)


tau = wline - 1  ##FIND A BETTER WAY!!!


for i in range(c.shape[1]):
    T = c[:, i].reshape((ms, ms))
    q = p[:, i]

    # print(T)
    # print(q)
    for k in range(ms):
        for l in range(ms):
            if q[k] > 0.0:
                L[k, l] = T[k, l] / q[k]
            else:
                L[k, l] = 0.0
                L[k, k] = 1.0

    # 	L = T / q[:,None]

    Ev, Evec = np.linalg.eig(np.transpose(L))
    idx = Ev.argsort()[::-1]
    Ev = Ev[idx]  # order eigenvectors by size3
    Evec = np.transpose(Evec)
    Evec = Evec[idx, :]
    tao = lagtime[i] * dT
    lt[0, i] = tao
    EV_arr[0, i] = tao  # first line of Evals is tao (for gnuplot)

    if i == tau:
        T_check = T
        np.savetxt(path0 + "single_particle/EV/ps_" + str(ident) + ".dat", q / sum(q))
        EV_comp[:, 0] = Ev.real
        EV_comp[:, 1] = Ev.imag

        Tplc = cpl[:, i].reshape((ms, ms))
        Tmic = cmi[:, i].reshape((ms, ms))
        Tpl = np.zeros((ms, ms))
        Tmi = np.zeros((ms, ms))
        for k in range(len(minpos)):
            target = list(range(minpos[k] - rancut, minpos[k] + rancut))
            for t in range(len(minpos)):
                origin = list(range(minpos[t] - rancut, minpos[t] + rancut))
                if origin[-1] >= ms or target[-1] >= ms:
                    print("origin/target for mfpt out of range!")
                else:
                    tau = wline * 1.0 - 1.0  # read from first line count...
                    mfpt_ar[t][k] = mfpt(L, target=target, origin=origin)

        MSM_ = MSM(L)

        for k in range(ms):
            for t in range(ms):
                # Probably not needed anymore
                # B = k
                # if (k < t):
                # 	B = [i+t for i,x in enumerate(r[k,t:]) if x == -1]
                # 	C = [i for i,x in enumerate(r[k,:k]) if x == -1]
                # 	B = B+C
                # 	if not B:
                # 		B = [i+k for i,x in enumerate(r[k,k:t+1]) if x == 1]
                # if (t < k):
                # 	B = [i for i,x in enumerate(r[k,:t+1]) if x == 1]
                # 	C = [i+k for i,x in enumerate(r[k,k:]) if x == 1]
                # 	B = B+C
                # 	if not B:
                # 		B = [i+t for i,x in enumerate(r[k,t:k]) if x == -1]
                # 		if ms in B:
                # 			B.remove(ms)

                tau = wline - 1.0
        # 				mfpt_arx[t][k] = mfpt(L, origin = k, target = t)

        length = 1000
        mfpt_ana = np.zeros((lvl * lvl, 1000))

        for k in range(lvl):
            for t in range(lvl):
                if k != t:
                    start = list(range(minpos[k] - rancut, minpos[k] + rancut + 1))
                    end = list(range(minpos[t] - rancut, minpos[t] + rancut + 1))
                    mfpt_ana[k * lvl + t, :] = first_passage(start, end, L, length)

        np.savetxt(
            path0 + "single_particle/MFPT/markovana_" + str(int(tau)) + "_" + str(ident) + ".dat",
            np.transpose(mfpt_ana),
        )
        # 		tra = MSM_.simulate(tra_len)
        # stuff = mfpt_trajectory_ms_cross(tra, blocks, minpos, cut, lvl, rancut, ms,r)

        # 		stuff = mfpt_trajectory_ms(tra, blocks, minpos, cut, lvl, rancut, ms)
        # 		MFPT_hist = np.zeros((9,1000))
        # 		MFPT_meas = stuff[0]
        # 		MFPT_err = stuff[1]
        # 		MFPT_steps = stuff[2]

        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/markov_'+str(int(tau))+"_"+str(ident)+'.dat',(mfpt_arx))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measma_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_meas))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measmaerr_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_err))

        # 		stuff = mfpt_trajectory_area(tra, blocks, minpos, cut, lvl, rancut, ms)

        MFPT_meas_reg = np.zeros((lvl, lvl))
        MFPT_err_reg = np.zeros((lvl, lvl))
        MFPT_steps_reg = np.zeros((lvl, lvl))
        MFPT_hist = np.zeros((lvl * lvl, 1000))
        # 		MFPT_meas_reg = stuff[0]
        # 		MFPT_err_reg = stuff[1]
        # 		MFPT_steps_reg = stuff[2]
        # 		MFPT_hist = stuff[3]
        # 		for j in range(lvl):
        # 			for k in range(lvl):
        # 				if (j !=k ):
        # 					MFPT_hist[j*lvl+k,:] = MFPT_hist[j*lvl+k,:] /sum(MFPT_hist[j*lvl+k,:])
        print(tau, ident)
        np.savetxt(path0 + "single_particle/MFPT/markov_area_" + str(int(tau)) + "_" + str(ident) + ".dat", (mfpt_ar))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measma_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_meas_reg))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measmaerr_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_err_reg))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/markovhist_'+str(int(tau))+"_"+str(ident)+'.dat',(np.transpose(MFPT_hist)))

        for k in range(ms):
            for l in range(ms):
                if q[k] > 0.0:
                    Tpl[k, l] = Tplc[k, l] / q[k]
                    Tmi[k, l] = Tmic[k, l] / q[k]
                else:
                    Tpl[k, l] = 0.0
                    Tmi[k, l] = 0.0
                    Tpl[k, k] = 1.0
                    Tmi[k, k] = 1.0

        np.savetxt(path0 + "single_particle/T/T_" + str(int(tau)) + "_" + str(ident) + ".dat", L)
        np.savetxt(path0 + "single_particle/T/Tpl_" + str(int(tau)) + "_" + str(ident) + ".dat", Tpl)
        np.savetxt(path0 + "single_particle/T/Tmi_" + str(int(tau)) + "_" + str(ident) + ".dat", Tmi)
        np.savetxt(path0 + "single_particle/EV/EV_com_" + str(int(tau)) + "_" + str(ident) + ".dat", EV_comp)
        for k in range(Enum):
            Evec_comp[:, 0] = Evec[k, :].real
            Evec_comp[:, 1] = Evec[k, :].imag
            for l in range(ms):
                if Evec[0, l] > 0:
                    Evecr_comp[l, 0] = (Evec[k, l] / Evec[0, l]).real
                    Evecr_comp[l, 1] = (Evec[k, l] / Evec[0, l]).imag
            # 				else:
            # 					Evecr_comp[l,0] = 0.
            # 					Evecr_comp[l,1] = 0.
            np.savetxt(path0 + "single_particle/EV/Evec" + str(k + 1) + "_com_" + str(ident) + ".dat", Evec_comp)
            np.savetxt(path0 + "single_particle/EV/Evecr" + str(k + 1) + "_com_" + str(ident) + ".dat", Evecr_comp)

    for j in range(1, ms):
        if Ev[j].real > 0 and not math.isclose(Ev[j].real, 1.0):  # How to deal with complex part?
            lt[j, i] = -tao / math.log(Ev[j].real)
        else:
            lt[j, i] = "NaN"
        EV_arr[j + 1, i] = Ev[j].real  # fill array with Ev's in order

    for k in range(Enum):
        Evec_arr[i, :, k] = Evec[k + 1].real  # Evec[0] cotains invariant distr.
        for l in range(ms):
            if Evec[0, l].real > 0:
                Evecr_arr[i, l, k] = Evec[k + 1, l].real / Evec[0, l].real  # right Evec
            else:
                Evecr_arr[i, l, k] = 0.0


if db_flag:  # do det bal test if wanted
    db_mat = np.zeros((ms, ms))
    av_violation = 0.0
    T_check = T_check.astype(float)

    T_check = T_check / np.sum(T_check[:, :])  # normalise

    for i in range(ms):
        for j in range(ms):
            db_mat[i, j] = abs(T_check[i, j] - T_check[j, i])
            av_violation += db_mat[i, j]
    av_violation /= (ms * ms - ms) / 2

    print("average violation of db: ", av_violation)

    np.savetxt(path0 + "single_particle/detbal/err_" + str(ident) + ".dat", db_mat)
    np.savetxt(path0 + "single_particle/detbal/flux_" + str(ident) + ".dat", T_check)


np.savetxt(path0 + "single_particle/lagtime/" + str(ident) + ".dat", np.transpose(lt))
np.savetxt(path0 + "single_particle/EV/EV_" + str(ident) + ".dat", np.transpose(EV_arr))


for k in range(Enum):
    np.savetxt(
        path0 + "single_particle/EV/Evec" + str(k + 1) + "_" + str(ident) + ".dat", np.transpose(Evec_arr[:, :, k])
    )

for k in range(Enum):
    np.savetxt(
        path0 + "single_particle/EV/Evecr" + str(k + 1) + "_" + str(ident) + ".dat", np.transpose(Evecr_arr[:, :, k])
    )
