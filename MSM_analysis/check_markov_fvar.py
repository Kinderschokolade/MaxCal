import argparse
import getpass
import math

import numpy as np
from define_flow import define_flow_direction, read_cut, read_minpos
from fct_markov import calc_av, calc_Sprod
from msmtools.analysis import mfpt
from scipy.optimize import fsolve
from scipy.special import erf

# from pyemma.msm import MSM
# from analyse_trajectory import mfpt_trajectory_ms_cross
# from analyse_trajectory import mfpt_trajectory_area
# from analyse_trajectory import mfpt_trajectory_ms
# from define_cross import define_cross
# from define_flow import define_flow_direction
from fpt_ana import first_passage


def escape_tr(L, q, c, relax_pos, lvl, ms):
    pi = np.zeros((lvl, ms))
    for j in range(lvl):
        for i in relax_pos[j]:
            pi[j, i] = q[i]
        pi[j] = pi[j] / sum(pi[j, :])  # starting vector

    p = np.zeros(2 * lvl)
    px = np.zeros((ms, lvl))

    for k in range(lvl):
        # transition matrix from  ms to A
        for l in relax_pos[k]:
            for m in relax_pos[k]:
                px[l, k] += L[l, m]
        # transition matrix from A to A
        for l in relax_pos[k]:
            p[k] += pi[k, l] * px[l, k]

        # err estmation from A to A
        p[lvl + k] = np.sqrt(
            k * (p[k] - p[k] * p[k]) / sum(sum(c[y * ms + z] for z in range(ms)) for y in relax_pos[k])
        )

    return p


def escape(pi, T, relax_pos, lvl, steps, lag):
    psum = np.zeros((lvl + 1, steps))
    ms = pi.size
    p = np.zeros((lvl, ms))
    p_old = np.zeros((lvl, ms))
    for k in range(lvl):
        for i in relax_pos[k]:
            psum[k + 1, 0] += pi[i]

    for j in range(lvl):
        for i in relax_pos[j]:
            p_old[j, i] = pi[i]
        p_old[j] = p_old[j] / sum(p_old[j, :])  # starting vector
        psum[j + 1, 0] = 1.0

    psum[0, 0] = 0

    for s in range(1, steps):
        # markov model estimation
        for k in range(lvl):
            for i in range(ms):
                for j in range(ms):
                    p[k, j] += p_old[k, i] * T[i, j]

        for k in range(lvl):
            for i in relax_pos[k]:
                psum[k + 1, s] += p[k, i]  # prob to stay in same basin
        psum[0, s] = s * lag
        p_old[:, :] = p[:, :]
        p[:, :] = 0.0

    return psum


def equ(c, ms, W):
    tup = [(sum([W[j, k] * c[k] * c[j] for k in range(ms)]) - 1.0) for j in range(ms)]
    tup = tuple(tup)
    return tup


def construct_MSM(T, r, Sprod, ms):
    gamma = 0.0  # for now
    eps = 1e-06
    W = np.zeros((ms, ms))
    for i in range(ms):
        for j in range(ms):
            if T[i, j] > eps and T[j, i] > eps:
                W[i, j] = T[i, j] * math.exp(
                    gamma / 2.0 * (r[i, j] + r[j, i]) + (Sprod[i, j] - np.log(T[i, j] / T[j, i])) / 2.0
                )
            else:
                W[i, j] = 1 * math.exp(gamma / 2.0 * (r[i, j] + r[j, i]) + (Sprod[i, j] - 1.0) / 2.0)

    c = np.ones(ms)

    copt = fsolve(
        equ,
        c,
        args=(
            ms,
            W,
        ),
        xtol=1e-8,
    )  # should converge for any starting point

    k = [[copt[j] * copt[i] * W[i, j] for j in range(ms)] for i in range(ms)]
    k = np.asarray(k)

    for i in range(ms):  # minimal correction from numerics, eta = 1.
        k[i, j] = k[i, j] / sum(k[i, :])

    return k


parser = argparse.ArgumentParser()
parser.add_argument("-o", help="ident")
parser.add_argument("-l", help="write T of specific line")
parser.add_argument("-n", help="number of (slowest) eigenvectors", default=3)
parser.add_argument("-db", help="set 1 for detbal check", default=0)
parser.add_argument("-gb", help="set 1 for global balance check", default=0)
parser.add_argument("-rc", help="rancut", default=1)
args = parser.parse_args()
ident = args.o
Enum = int(args.n)
tau = int(args.l)
db_flag = int(args.db)
global_flag = int(args.gb)
rancut = int(args.rc)

user = getpass.getuser()

if user == "mbause":
    path0 = "/u/mbause/data/"
    path1 = "/u/mbause/data/"
else:
    path0 = "/data/isilon/bause/"
    path1 = "/data/pckr194/bause/"


string = " "
ms = 0
dT = 0
for line in open(path0 + "single_particle/param_" + str(ident) + ".dat", "r"):
    cont = line.split()
    if cont[0] == "dT":
        dT = float(cont[1])
    if cont[0] == "microstates":
        ms = int(cont[1])
    if cont[0] == "T":
        Temperature = float(cont[1])
    if cont[0] == "fmax":
        fmax = float(cont[1])
    if cont[0] == "fsig":
        fsig = float(cont[1])
    if cont[0] == "fx":
        fx = float(cont[1])


c_s = np.zeros(ms * ms)
p_s = np.zeros(ms)

F = open(path0 + "single_particle/count_" + str(ident) + ".dat", "r")
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime = [int(i) for i in lagtime]
print(lagtime)

cols = len(lagtime)

data = np.loadtxt(path0 + "single_particle/count_" + str(ident) + ".dat", dtype="i", usecols=(range(1, cols + 1)))


###########
potential = np.loadtxt(path0 + "single_particle/potential_ms_" + str(ident) + ".dat")
###############

ms, minpos = read_minpos(path0, ident)
ms, cut = read_cut(path0, ident)
F = define_flow_direction(ms, cut)

c = data[0 : ms * ms]
p = data[ms * ms : ms * ms + ms]

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

force = np.zeros(ms)
for i in range(ms):
    pos = (i + 1.0 / 2.0) / ms - fx
    force[i] = fmax * i

lvl = len(minpos)

rc_bas = 4
relax_pos = [[0] * (rc_bas * 2 + 1)] * lvl
for i in range(lvl):
    relax_pos[i] = list(range(minpos[i] - rc_bas, minpos[i] + rc_bas + 1))
    for j in range(rc_bas * 2 + 1):
        if relax_pos[i][j] < 0:
            relax_pos[i][j] = ms - relax_pos[i][j]
        if relax_pos[i][j] > (ms - 1):
            relax_pos[i][j] = relax_pos[i][j] - ms


mfpt_ar = np.zeros((lvl, lvl))
mfpt_arx = np.zeros((ms, ms))
mfpt_meas = np.zeros((ms, ms))
mfpt_err = np.zeros((ms, ms))
mfpt_steps = np.zeros((ms, ms))


blocks = 10
tra_len = 1000000

# r = flow_direction(cpl[:,wline].reshape((ms,ms)), cmi[:,wline].reshape((ms,ms)), ms)
Sp = np.zeros((ms, ms))
for k in range(ms):
    for l in range(ms):
        cutoff = 0.0001
        if np.abs(k - l) < (ms / 2):  # why do I change sign here (compared to reweighting)
            if l > k:
                Sp[k, l] = (np.sqrt(math.pi / 2.0) * (erf()) + potential[k, 1] - potential[l, 1]) / Temperature
            else:
                Sp[k, l] = (+sum(force[l:k]) + potential[k, 1] - potential[l, 1]) / Temperature
        elif l > k:
            Sp[k, l] = (sum(force[l:ms]) + sum(force[0:k]) + potential[k, 1] - potential[l, 1]) / Temperature
        else:
            Sp[k, l] = (-sum(force[k:ms]) - sum(force[0:l]) + potential[k, 1] - potential[l, 1]) / Temperature
# 		print(k,l, Sp[k,l])


tao_MFPT = np.zeros((c.shape[1], 10))
tao_J = np.zeros((c.shape[1], 2))
tao_Sprod = np.zeros((c.shape[1], 2))

steps = 100
p_ck = np.zeros((lvl + 1, steps))
p_tr = np.zeros((2 * lvl + 1, c.shape[1]))

# escape time trajectory based (ck-test)

for i in range(cols):
    T = np.transpose(c[:, i].reshape((ms, ms)))
    listdel = []
    q = p[:, i]
    for k in range(ms):
        if q[k] < 0.1:
            listdel.append(k)

    for k in range(ms):
        for l in range(ms):
            if q[k] > 0.0:
                L[k, l] = T[k, l] / sum(T[k, :])  # / q[k]
            else:
                L[k, l] = 0.0
                L[k, k] = 1.0

    Lred = np.delete(L, listdel, axis=0)
    Lred = np.delete(Lred, listdel, axis=1)
    # print(lagtime[i])
    # print(Lred)

    Ev, Evec = np.linalg.eig(np.transpose(Lred))
    idx = Ev.argsort()[::-1]
    Ev = Ev[idx]  # order eigenvectors by size
    Evec = np.transpose(Evec)
    Evec = Evec[idx, :]
    p_stat = np.real(Evec[0, :] / sum(Evec[0, :]))
    tao = lagtime[i]
    lt[0, i] = tao
    EV_arr[0, i] = tao  # first line of Evals is tao (for gnuplot)
    tao_MFPT[i][0] = tao

    p_tr[0, i] = tao
    p_tr[1:, i] = escape_tr(L, q, c[:, i], relax_pos, lvl, ms)

    # 	for k in range(len(minpos)):
    # 		target = list(range(minpos[k] - rancut , minpos[k] + rancut +1))
    # 		for t in range(len(minpos)):
    # 			origin = list(range(minpos[t] - rancut , minpos[t] + rancut +1))
    # 			if(origin[-1] >= ms or target[-1] >= ms):
    # 				print ("origin/target for mfpt out of range!")
    # 			else:
    # 				tao_MFPT[i][1+t*3+k] = mfpt(L,target = target, origin = origin)

    tao_J[i, 0] = tao
    tao_Sprod[i, 0] = tao
    tao_J[i, 1] = 0.0  # calc_av(F, L, p_stat)  # error because p_stat is of antoher size
    tao_Sprod[i, 1] = calc_Sprod(L)

    if lagtime[i] == tau:
        L_check = L
        #################################
        # L = construct_MSM(L,F,Sp,ms)
        #################################
        T_check = T

        np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/EV/ps_" + str(ident) + ".dat", q / sum(q))
        p_check = q / sum(q)
        # EV_comp[:,0] =  Ev.real # size problem
        # EV_comp[:,1] =  Ev.imag

        Tc = np.transpose(c[:, i].reshape((ms, ms)))
        # 		np.savetxt("count.dat",Tplc+Tmic)
        T = np.zeros((ms, ms))
        # mean first passage time
        for k in range(len(minpos)):
            target = list(range(minpos[k] - rancut, minpos[k] + rancut + 1))
            for t in range(len(minpos)):
                origin = list(range(minpos[t] - rancut, minpos[t] + rancut + 1))
                if origin[-1] >= ms or target[-1] >= ms:
                    print("origin/target for mfpt out of range!")
        # 				else:
        # 					tau = wline * 1. - 1.  # read from first line count...
        # 					mfpt_ar[t][k] = mfpt(L,target = target, origin = origin)

        # escapetime for MSM chapman kolmogorov test
        p_ck = escape(p_check, L, relax_pos, lvl, steps, tau)

        # MSM_ = MSM(L)

        length = 1000
        fpt_ana = np.zeros((lvl * lvl, length))

        # First passage Time
        for k in range(lvl):
            for t in range(lvl):
                if k != t:
                    start = list(range(minpos[k] - rancut, minpos[k] + rancut + 1))
                    end = list(range(minpos[t] - rancut, minpos[t] + rancut + 1))
                    fpt_ana[k * lvl + t, :] = first_passage(start, end, L, length)
                    # print(start,end)

        np.savetxt(
            path0 + "single_particle/" + str(int(tau)) + "/MFPT/markovana_" + str(ident) + ".dat", np.transpose(fpt_ana)
        )
        np.savetxt(
            path0 + "single_particle/" + str(int(tau)) + "/MFPT/cktest_" + str(ident) + ".dat", np.transpose(p_ck)
        )

        ratio = np.zeros((ms, ms))
        ratio_check = np.zeros((ms, ms))
        for k in range(ms):
            for l in range(ms):
                cutoff = 0.0001
                if L[k, l] > cutoff and L[l, k] > cutoff:
                    ratio[k, l] = np.log(L[k, l] / L[l, k])
                if L_check[k, l] > cutoff and L_check[l, k] > cutoff:
                    ratio_check[k, l] = np.log(L_check[k, l] / L_check[l, k])

        print("/data/isilon/bause/single_particle/" + str(int(tau)) + "/T/ratio_" + str(ident) + ".dat")

        np.savetxt("/data/isilon/bause/single_particle/" + str(int(tau)) + "/T/ratio_" + str(ident) + ".dat", ratio)
        np.savetxt(
            "/data/isilon/bause/single_particle/" + str(int(tau)) + "/T/ratio_check_" + str(ident) + ".dat", ratio_check
        )
        np.savetxt("/data/isilon/bause/single_particle/" + str(int(tau)) + "/T/ratio_th_" + str(ident) + ".dat", Sp)

        # test FPT by comparing to trajectory prododuced by MSM
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
        # np.savetxt(path0+'single_particle/'+str(int(tau))+'/MFPT/markov_area_'+str(ident)+'.dat',(mfpt_ar))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measma_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_meas_reg))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/measmaerr_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_err_reg))
        # 		np.savetxt('/data/isilon/bause/single_particle/MFPT/markovhist_'+str(int(tau))+"_"+str(ident)+'.dat',(np.transpose(MFPT_hist)))

        for k in range(ms):
            for l in range(ms):
                if q[k] > 0.0:
                    T[k, l] = Tc[k, l] / q[k]
                else:
                    T[k, l] = 0.0
                    T[k, k] = 1.0

        np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/T/T_" + str(ident) + ".dat", L)
        np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/T/Tcl_" + str(ident) + ".dat", L_check)

        np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/EV/EV_com_" + str(ident) + ".dat", EV_comp)

        Omega = np.zeros((ms, ms))
        for k in range(ms):
            for t in range(ms):
                Omega[k, t] = np.sqrt(
                    L[k, t] * L[k, t]
                )  # * np.exp((potential[k,1] + potential[t,1] )/(-2*Temperature))

        np.savetxt("/data/isilon/bause/single_particle/" + str(int(tau)) + "/T/Omega_" + str(ident) + ".dat", Omega)
        # size problem
        # for k in range(Enum):
    # 			Evec_comp[:,0] = Evec[k,:].real
    # 			Evec_comp[:,1] = Evec[k,:].imag
    # 			for l in range(ms):
    # 				if (Evec[0,l] > 0):
    # 					Evecr_comp[l,0] = (Evec[k,l]/ Evec[0,l]).real
    # 					Evecr_comp[l,1] = (Evec[k,l]/ Evec[0,l]).imag
    # 				else:
    # 					Evecr_comp[l,0] = 0.
    # 					Evecr_comp[l,1] = 0.
    # 			np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/Evec'+str(k+1)+'_com_'+str(ident)+'.dat',Evec_comp)
    # 			np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/Evecr'+str(k+1)+'_com_'+str(ident)+'.dat',Evecr_comp)

    msred = Lred.shape[0]

    for j in range(1, msred):
        if Ev[j].real > 0 and not math.isclose(Ev[j].real, 1.0):  # How to deal with complex part?
            lt[j, i] = -tao / math.log(Ev[j].real)
        else:
            lt[j, i] = "NaN"
        EV_arr[j + 1, i] = Ev[j].real  # fill array with Ev's in order
# size problem
# 	for k in range(Enum):
# 		Evec_arr[i,:,k] = Evec[k+1].real # Evec[0] cotains invariant distr.
##		for l in range(msred):
##			if (Evec[0,l].real > 0):
# 				Evecr_arr[i,l,k] = Evec[k+1,l].real / Evec[0,l].real #right Evec
# 			else:
# 				Evecr_arr[i,l,k] = 0.


np.savetxt(path0 + "single_particle/lagtime/J_" + ident + ".dat", tao_J)
np.savetxt(path0 + "single_particle/lagtime/ck_" + ident + ".dat", np.transpose(p_tr))
np.savetxt(path0 + "single_particle/lagtime/Sprod_" + ident + ".dat", tao_Sprod)
np.savetxt(path0 + "single_particle/lagtime/MFPT_" + ident + ".dat", tao_MFPT)

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

    np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/detbal/err_" + str(ident) + ".dat", db_mat)
    np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/detbal/flux_" + str(ident) + ".dat", T_check)


if global_flag:  # do global bal test if wanted
    err = np.zeros(ms)
    glob = np.zeros(ms)

    for i in range(ms):
        err[i] = abs(sum(p_check[j] * L_check[j, i] for j in range(ms)) - p_check[i])
        glob[i] = abs(sum(p_check[j] * L_check[j, i] for j in range(ms)))

    av_violation = sum(err) / ms

    p_check = np.transpose([p_check])
    glob = np.transpose([glob])

    out = np.concatenate((glob, p_check), axis=1)

    np.savetxt(path0 + "single_particle/" + str(int(tau)) + "/detbal/glob_" + str(ident) + ".dat", out)
    print("average violation of global balance: ", av_violation)


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
