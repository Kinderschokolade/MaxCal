import numpy as np


def stat_distribution(T):
    Ev, Evec = np.linalg.eig(np.transpose(T))
    idx = Ev.argsort()[::-1]
    Ev = Ev[idx]  # order eigenvectors by size
    Evec = np.transpose(Evec)
    Evec = Evec[idx, :]
    p = np.real(Evec[0] / sum(Evec[0]))
    return p


def first_passage(start, end, T, length):
    # start and end are single ms for now
    n = T.shape[0]
    T_a = np.zeros((n, n))
    T_a[:, :] = T

    if isinstance(end, list):
        N_end = len(end)
    else:
        N_end = 1
        end = [end]

    if isinstance(start, list):
        N_start = len(start)
    else:
        N_start = 1
        start = [start]

    # calc by choosing absorbing points in end and starting with probability 1 at starting point
    fp = np.zeros((N_start, length))
    for i in range(N_end):
        T_a[end[i], :] = 0
        T_a[end[i], end[i]] = 1
        # final points are now absorbing

    J = np.zeros((N_start, n))
    for j in range(N_start):  # perform for each starting point
        pi = np.zeros((length + 1, n))
        prob = np.zeros(length + 1)
        pi[0, start[j]] = 1
        for t in range(length):
            for i in range(N_end):
                prob[t] += pi[t, end[i]]  # all prob of being in state end are collected in prob
            pi[t + 1, :] = pi[t, :].dot(T_a)  # update to new conf
            # iterate Markov process

        fp[j, 1:length] = [prob[i] - prob[i - 1] for i in range(1, length)]
        fp[j, 0] = prob[0]
    # 		for a in range(n):
    # 			if (a not in start and a not in end):
    # 				J[j,a] += sum(pi[:,a])  ##needs weights correction later
    # use another time!!!!!!!!!!!!!

    # calc stationary distribution
    Ev, Evec = np.linalg.eig(np.transpose(T))
    idx = Ev.argsort()[::-1]
    Ev = Ev[idx]  # order eigenvectors by size
    Evec = np.transpose(Evec)
    Evec = Evec[idx, :]
    p = np.real(Evec[0] / sum(Evec[0]))
    fp_full = np.zeros(length)

    # weight by prob that TRAJECTORY starts in start
    weight = np.zeros(n)
    for j in start:
        for k in range(n):
            if k not in start:
                weight[j] += p[k] * T[k, j]

    weight = weight / sum(weight)  # normalise

    # 	J_full = np.zeros(n)
    for j in range(N_start):
        fp_full += weight[start[j]] * fp[j]
    # 		J_full += weight[start[j]] * J[j]

    fp_full = fp_full / sum(fp_full)  # corrects for very small numercial error
    # 	J_full= J_full / sum(J_full)

    # 	return fp_full, J_full
    return fp_full
