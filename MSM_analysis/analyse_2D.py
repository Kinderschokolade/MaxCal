import math

import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy import linalg
from scipy.optimize import fmin, fsolve

############################################

# BASICS

############################################


##### returns data for a data index #####
def get_data(out, data_index):  # out = 1,2,3 gives count_matrix, potential or parameters respectively
    path = "/data/isilon/bause/single_particle/"
    if out == 1:
        return (path + "count_{}.dat").format(data_index)
    if out == 2:
        data_index = str(data_index)[0:2] + "0" + str(data_index)[3]
        return (path + "potential_{}.dat").format(data_index)
    if out == 3:
        data_index = str(data_index)[0:2] + "0" + str(data_index)[3]
        return (path + "param_{}.dat").format(data_index)


##### get the shape of the MSM model, i.e. [30,10] for 2D model #####
def get_shape(data_index):
    data = [get_data(3, data_index)]
    with open(data[0]) as dat:  # get params from param file
        for line in dat:
            if line.split("	")[0] == "microstates ":
                num_states = line.split("	")[1].split(" ")[::2]
                num_states = [int(i) for i in num_states]
    return num_states


##### get the external force strength of the simulation, i.e. [30,10] for 2D model #####
def get_force(data_index):
    data = [get_data(3, data_index)]
    with open(data[0]) as dat:  # get params from param file
        for line in dat:
            if line.split("	")[0] == "extf0":
                f = float(line.split("	")[2])
    return f


##### get the transition matrix for a certain data_index and lagtime #####
def get_transition_matrix(data_index, lagtime):
    shape = get_shape(data_index)
    data = [get_data(1, data_index)]
    count_matrix = np.zeros((shape[0] * shape[1], shape[0] * shape[1]))
    transition_matrix = np.zeros((shape[0] * shape[1], shape[0] * shape[1]))
    with open(data[0]) as dat:  # put all data in list p
        p = []
        for line in dat:
            p.append(float(line.split("	")[lagtime]))
    for i in range(0, shape[0] * shape[1]):  # assigns data to the count matrix entries
        for j in range(0, shape[0] * shape[1]):
            count_matrix[i][j] = p[shape[0] * shape[1] * i + j + 1]

    for i in range(0, shape[0] * shape[1]):  # get transition matrix
        for j in range(0, shape[0] * shape[1]):
            if p[i + shape[0] * shape[1] * shape[0] * shape[1] + 1] != 0:
                transition_matrix[i][j] = (
                    count_matrix[j][i] / p[i + shape[0] * shape[1] * shape[0] * shape[1] + 1]
                )  # count matrix devided by the row sum
            else:
                transition_matrix[i][j] = 0
                transition_matrix[i][i] = 1
    return transition_matrix


##### prints parameter file #####
def get_params(data_index):
    data = [get_data(3, data_index)]
    with open(data[0]) as dat:
        for line in dat:
            print(line)


##### returns the stationary distribution of a given transition matrix #####
def get_stat_distribution(transition_matrix):
    w, vec = np.linalg.eig(transition_matrix.T)  # transpose matrix to get the left eigenvectors
    idx = w.argsort()[::-1]
    w = w[idx]  # order eigenvectors by size
    vec = np.transpose(vec)
    vec = vec[idx, :]
    p = np.real(vec[0] / sum(vec[0]))
    return p


##### plots the 2D potential #####
def plot_potential(data_index, figax=False):  # if figax = True returns fig,ax
    data = [get_data(2, data_index)]

    shape = get_shape(data_index)
    N = shape[0] * shape[1]

    x_data = np.fromfile(data[0], sep="	")[0::N]
    y_data = np.fromfile(data[0], sep="	")[1:N:3]
    pot_data = np.fromfile(data[0], sep="	")[2::3]

    X, Y = np.meshgrid(x_data, y_data)
    pot_data = pot_data.reshape(X.shape).T

    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(111, projection="3d")
    p = ax.plot_wireframe(X, Y, pot_data)
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel("$y$", fontsize=12)
    ax.set_zlabel("$V(x,y)$", fontsize=12)
    plt.show()

    if figax == True:
        return fig, ax


############################################

# MSM ANALYSIS

############################################


##### returns eigenvalues and both eigenvectors with seperated real and imaginary part #####
def get_eigen(data_index, lagtime):
    transition_matrix = get_transition_matrix(data_index, lagtime)
    w, vl, vr = linalg.eig(transition_matrix, left=True, right=True)

    eigen_lst = []
    for i in range(len(w)):
        eigen_lst.append((np.real(w[i]), np.real(vl[:, i]), np.imag(vl[:, i]), np.real(vr[:, i]), np.imag(vr[:, i])))
    eigen_lst = sorted(
        eigen_lst, key=lambda tup: tup[0], reverse=True
    )  # sort everything by the size of the eigenvalues

    w = [i[0] for i in eigen_lst]  # eigenvalues

    vl = ([i[1] for i in eigen_lst], [i[2] for i in eigen_lst])  # eigenvectors with seperated real and imag part
    if vl[0][0][0] < 0:  # the first eigenvectors are now positive (just for visualization)
        vl[0][0] = -vl[0][0]
    vr = ([i[3] for i in eigen_lst], [i[4] for i in eigen_lst])
    if vr[0][0][0] < 0:
        vr[0][0] = -vr[0][0]

    return w, vl, vr


##### plots eigenvalues #####
def plot_eigenvalues(data_index, lagtime, figax=False):  # if figax = True returns fig,ax
    shape = get_shape(data_index)
    w, vl, vr = get_eigen(data_index, lagtime)
    x_axis = np.arange(1, shape[0] * shape[1] + 1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x_axis, w, linestyle="None", marker="x", markersize=8)
    plt.show()
    if figax == True:
        return fig, ax


##### plots eigenvectors with seperated real and imaginary part #####
def plot_eigenvectors(
    data_index, lagtime, ew=1, left=True, right=False, imag=False
):  # choose which eigenvectors (ew=?), right\left shall be plotted and if the imaginary parts shall be plotted aswell
    shape = get_shape(data_index)
    w, vl, vr = get_eigen(data_index, lagtime)

    if left == True:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        x = np.arange(1, shape[0] + 1)
        y = np.arange(1, shape[1] + 1)
        X, Y = np.meshgrid(x, y)
        ax.plot_surface(X, Y, vl[0][ew].reshape(X.shape))
        plt.show()
        if imag == True:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection="3d")
            x = np.arange(1, shape[0] + 1)
            y = np.arange(1, shape[1] + 1)
            X, Y = np.meshgrid(x, y)
            ax.plot_surface(X, Y, vl[1][ew].reshape(X.shape))
            plt.show()

    if right == True:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        x = np.arange(1, shape[0] + 1)
        y = np.arange(1, shape[1] + 1)
        X, Y = np.meshgrid(x, y)
        ax.plot_surface(X, Y, vr[0][ew].reshape(X.shape))
        plt.show()
        if imag == True:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection="3d")
            x = np.arange(1, shape[0] + 1)
            y = np.arange(1, shape[1] + 1)
            X, Y = np.meshgrid(x, y)
            ax.plot_surface(X, Y, vr[1][ew].reshape(X.shape))
            plt.show()


##### returns lagtime-timescale plot for different values #####
def plot_lagtime(
    data_index, values, max_time=999, step=10, figax=False
):  # values has to be a list of the different eigenvalues which you want the plot for, if figax = True returns fig,ax
    fig, ax = plt.subplots(figsize=(10, 6))

    tau_lst = []
    for i in range(0, max_time, step):  # lagtime list
        tau = i + 1
        tau_lst.append(tau)

    for value in values:
        ti_lst = []

        for i in range(0, max_time, step):  # calculate timescale
            tau = (i + 1) * 1e-5
            w, vl, vr = get_eigen(
                data_index, i + 1
            )  ###TO DO not very efficient to calculate vectors aswell since they are not needed
            ti = -tau / np.log(w[value])
            ti_lst.append(ti)

        ax.plot(tau_lst, ti_lst, label="t{}".format(value), marker=".", markersize=2, linestyle="")

    ax.legend(loc="best", fontsize=12)
    ax.set_xlabel("tau", fontsize=12)
    ax.set_ylabel("$t_i$", fontsize=12)
    plt.show()

    if figax == True:
        return fig, ax


##### performs Chapman Kolmogorov test, trajectory needed #####
def ck_test(data_index, lag, set_A, k_range, plot=False):
    # set_A has to be a list of states which resemble a metastable macrostate
    # k_range the range of k values we want to look at

    shape = get_shape(data_index)

    # get trajectory in terms of microstates (path has to be adjusted)
    for i in range(20):
        ident = str(data_index) + "_{}".format(i)

        source = open("/home/../data/isilon/bause/single_particle/trajectory_" + str(ident) + ".dat", "r")

        num = 10000000
        for i in range(num):
            stuff = source.readline()

            if len(stuff) > 1:
                x = int(float(stuff.split()[1]) * shape[1])
                y = int(float(stuff.split()[2]) * shape[1])
                traj.append(y * shape[0] + x)

    T = get_transition_matrix(data_index, lag)
    pi = get_stat_distribution(T)

    # get the count matrix from the trajectory for different windows k
    def get_count(traj, k):
        C = np.zeros((shape[0] * shape[1], shape[0] * shape[1]))
        traj = traj[:: k * lag]
        for i in range(len(traj) - 1):
            C[traj[i]][traj[i + 1]] += 1
        return C

    # probability to still be in set_A after k steps according to the trajectory p_MD
    def w_md(i):
        if i in set_A:
            return pi[i] / sum([pi[j] for j in set_A])
        else:
            return 0

    def pmd(i, k):
        C = get_count(traj, k)
        return sum([C[i][j] for j in set_A]) / sum([C[i][j] for j in range(shape[0] * shape[1])])

    def p_MD(k):
        return sum([w_md(i) * pmd(i, k) for i in set_A])

    def err_MD(k):
        C = get_count(traj, k)
        return np.sqrt(
            k * (p_MD(k) - p_MD(k) ** 2) / sum([sum([C[i][j] for j in range(shape[0] * shape[1])]) for i in set_A])
        )

    # probability to still be in set_A after k steps according to the MSM p_MSM
    w = np.zeros(shape[0] * shape[1])
    for i in range(shape[0] * shape[1]):
        w[i] = w_md(i)

    def p_MSM(k):
        Tk = np.eye(shape[0] * shape[1])
        for i in range(k):
            Tk = np.dot(Tk, T)
        return sum([np.dot(w, Tk)[i] for i in set_A])

    # gives list of probabilities for MSM, traj and traj_err for different k values
    MSM_ck = [p_MSM(k) for k in k_range]
    MD_ck = [p_MD(k) for k in k_range]
    MD_ck_err = [err_MD(k) for k in k_range]

    # plots ck_test
    if plot == True:
        fig, ax = plt.subplots(dpi=150)
        ax.plot(k_range, MSM_ck, linestyle="", marker=".", label="Markov state model")
        ax.errorbar(k_range, MD_ck, MD_ck_err, linestyle="", capsize=3, marker="D", markersize=3, label="trajectory")
        ax.set_xlabel("k [tau]")
        ax.set_ylabel("p")
        ax.legend(loc="best")
        plt.show()

    return MSM_ck, MD_ck, MD_ck_err


##### show histogramm of count matrix #####
def show_matrix(data_index, lagtime, potential=True):  # if potential = True also plot the potential
    shape = get_shape(data_index)
    N = shape[0] * shape[1]
    data = [get_data(1, data_index), get_data(2, data_index)]
    with open(data[0]) as dat:  # plot histogramm matrix
        p = []
        for line in dat:
            p.append((line.split("	")[lagtime]))
        H = np.array([float(i) for i in p[(shape[0] * shape[1]) ** 2 + 1 :]])
    H = H.reshape(shape[1], shape[0])
    fig, ax = plt.subplots(figsize=(20, 20))
    p = ax.imshow(H, interpolation="nearest", origin="lower", extent=[0, 3, 0, 1])
    fig.colorbar(p, shrink=0.3)
    plt.show()

    if potential == True:  # plot potential
        x_data = np.fromfile(data[1], sep="	")[0::N]
        y_data = np.fromfile(data[1], sep="	")[1:N:3]
        pot_data = np.fromfile(data[1], sep="	")[2::3]
        X, Y = np.meshgrid(x_data, y_data)
        pot_data = pot_data.reshape(X.shape)

        q = ax.contour(X, Y, pot_data.T, cmap="winter")
        fig.colorbar(q, shrink=0.3)
        plt.show()


##### visulaize the transition matrix in terms of a starting state 'ref_state' #####
def visualize_Tij(data_index, lagtime, ref_state):

    T = get_transition_matrix(data_index, lagtime)
    shape = get_shape(data_index)

    for i in range(len(T[ref_state])):
        if T[ref_state][i] == 0:
            T[ref_state][i] = None
    fig, ax = plt.subplots(dpi=200)
    p = ax.imshow(T[ref_state].reshape(shape[1], shape[0]), interpolation="nearest", origin="lower", cmap="Reds")
    fig.colorbar(p, tick_params(labelsize=6), shrink=0.3)
    ax.plot(10, 4, marker="x", markersize=8, linestyle="", label="reference point")
    ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("No. of micro-state on x-axis", fontsize=8)
    ax.set_ylabel("No. of micro-state on y-axis", fontsize=8)
    ax.tick_params(axis="both", labelsize=8)


############################################

# REWEIGHTING

############################################


##### defines the flux matrix #####
def get_flux_matrix(shape):
    flux = np.zeros([shape[0] * shape[1], shape[0] * shape[1]])
    for line in range(len(flux)):
        minus = np.mod(line, shape[0])
        prior = [0] + list(np.arange(1, shape[0] / 2)) + [0] + list(np.arange(-shape[0] / 2 + 1, 0))
        seq = prior[shape[0] - minus :] + prior[: shape[0] - minus]
        flux[line] = np.array(seq * shape[1]) * 3 / shape[0]  ### distance metric  ### 3 = x_length
    return flux


##### get the average flux from a simulation transition matrix #####
def average_flux(data_index, lagtime):
    shape = get_shape(data_index)
    flux = get_flux_matrix(shape)
    trans = get_transition_matrix(data_index, lagtime)
    p = get_stat_distribution(trans)

    matrix = np.zeros([shape[0] * shape[1], shape[0] * shape[1]])
    for a in range(len(matrix)):
        for b in range(len(matrix[a])):
            matrix[a][b] = p[a] * trans[a][b] * flux[a][b] / lagtime / 1e-5  ## scaling flux with lagtime

    return sum(sum(matrix))


##### help funtion which returns the functions to determine the reweighted matrix and its flux only depending on gamma #####
def rew_functions(data_index, lagtime):  ##saves computing time if individual functions dont need to acces shape and q
    shape = get_shape(data_index)
    flux = get_flux_matrix(shape)
    q = get_transition_matrix(data_index, lagtime)

    def get_rew_matrix(y):

        ## matrix W
        W = q * np.exp(y / lagtime * flux)

        ## get first eigenvalue and vector of W
        w, v = linalg.eig(W)
        idx = w.argsort()[::-1]
        w = w[idx]
        v = np.transpose(v)
        v = v[idx, :]
        w, v = np.real(w[0]), np.real(v[0])

        factor = np.zeros([shape[0] * shape[1], shape[0] * shape[1]])
        for a in range(len(factor)):
            for b in range(len(factor[a])):
                factor[a][b] = v[b] / (v[a] * w)

        return factor * W

    def get_rew_flux(y):

        matrix = np.zeros([shape[0] * shape[1], shape[0] * shape[1]])
        k_matrix = get_rew_matrix(y)
        pi_vec = get_stat_distribution(k_matrix)

        for a in range(len(matrix)):
            for b in range(len(matrix[a])):
                matrix[a][b] = pi_vec[a] * k_matrix[a][b] * flux[a][b] / lagtime / 1e-5  ## scaling flux with lagtime

        return sum(sum(matrix))

    return get_rew_matrix, get_rew_flux


##### reweighting from a reference model to a target system #####
def reweighting_old(ref_data, target_data, lagtime):
    target_flux = average_flux(target_data, lagtime)
    rew_matrix, rew_flux = rew_functions(ref_data, lagtime)

    gamma = fmin(lambda y: (rew_flux(y) - target_flux) ** 2, 0, disp=False)

    return rew_matrix(gamma)  # returns reweighted transition matrix


#### NEW REWEIGHTING


##### calculates the local entropy production from the model as well as the theoretical value from state 'start' to state 'end' #####
def entroprod_matrix(data_index):

    potentialeval, defdomain, num_states = Gauss_potential(data_index)  # get the potential function
    x_min, x_max, y_min, y_max = defdomain

    N = num_states[0] * num_states[1]

    # defines microstates, 'states' describes a list of all microstates defined by their center
    states_h = [
        np.arange(x_min + 1 / (2 * num_states[0]) * (x_max - x_min), x_max, (x_max - x_min) / num_states[0]),
        np.arange(y_min + 1 / (2 * num_states[1]) * (y_max - y_min), y_max, (y_max - y_min) / num_states[1]),
    ]
    states_h = [[[j, i] for j in states_h[0]] for i in states_h[1]]
    states = []
    for i in states_h:
        for j in i:
            states.append(j)

    f = get_force(data_index)  # external force strength (only works for a constant force)
    T = 5  # need to be careful here if you change the temperature T

    S = np.zeros((N, N))

    for start in range(N):
        for end in range(N):

            # calculate distance between the two states, accounts for the direction of the transition in question
            val = states[start][0] - states[end][0]
            if val > x_max / 2:
                dist = val - x_max
            elif val < -x_max / 2:
                dist = x_max + val
            elif val == x_max / 2 or val == -x_max / 2:
                dist = 0
            else:
                dist = val
            S[start][end] = (
                -f * dist * 0.97 + potentialeval(states[start]) - potentialeval(states[end])
            ) / T  # 0.945 is a manual factor

    return S


def reweighting(ref_data, target_data, lagtime):

    q = get_transition_matrix(ref_data, lagtime)
    N = len(q)

    S_ref = entroprod_matrix(ref_data)
    S_tar = entroprod_matrix(target_data)

    W = q * np.exp(0.5 * (S_tar - S_ref))

    Phi = np.ones(N)

    Phi_opt = fsolve(lambda Phi: [np.dot(W, Phi)[i] * Phi[i] - 1.0 for i in range(N)], Phi)

    p = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            p[i][j] = W[i][j] * Phi_opt[i] * Phi_opt[j]

    for i in range(N):
        for j in range(N):
            p[i, j] = p[i, j] / sum(p[i, :])

    return p


############################################

# FIRST PASSAGE TIME DISTRIBUTIONS

############################################


##### defines the Gaussian potential #####
def Gauss_potential(data_index):

    data = [get_data(3, data_index)]
    with open(data[0]) as dat:  # get params from param file
        for line in dat:
            if line.split("	")[0] == "microstates ":
                num_states = line.split("	")[1].split(" ")[::2]
                num_states = [int(i) for i in num_states]
            if line.split("	")[0] == "sigma1 ":
                sigma1 = float(line.split("	")[1])
            if line.split("	")[0] == "sigma1":
                sigma1 = float(line.split("	")[2])
            if line.split("	")[0] == "sigma2 ":
                sigma2 = float(line.split("	")[1])
            if line.split("	")[0] == "sigma2":
                sigma2 = float(line.split("	")[2])
            if line.split("	")[0] == "U0":
                U0 = float(line.split("	")[2])
            if line.split("	")[0] == "U1":
                U1 = float(line.split("	")[2])
            if line.split("	")[0] == "U2":
                U2 = float(line.split("	")[2])

    defdomain = (0, 3, 0, 1)  # domain of definiton for the potential (has to be adjusted manually)

    def potentialeval(pos):  # 3 Gaussian wells at different positions
        Gauss0 = -U0 * np.exp(-((pos[0] - 0.5) ** 2 / 2 / sigma1**2 + (pos[1] - 0.5) ** 2 / 2 / sigma2**2))
        Gauss1 = -U1 * np.exp(-((pos[0] - 1.5) ** 2 / 2 / sigma1**2 + (pos[1] - 0.5) ** 2 / 2 / sigma2**2))
        Gauss2 = -U2 * np.exp(-((pos[0] - 2.5) ** 2 / 2 / sigma1**2 + (pos[1] - 0.5) ** 2 / 2 / sigma2**2))
        return Gauss0 + Gauss1 + Gauss2

    return potentialeval, defdomain, num_states


##### assignes microstates to the different potential wells #####
def state_assignement(potential, data_index, threshhold=0.1, disp=True):
    # potential has to be a 2D funtion in the form of the Gauss_potential - function,
    # the threshhold defines until which energy level in the well microstates will get assigned to it
    # if disp = True: prints the assigned states and determined potential wells

    potentialeval, defdomain, num_states = potential(data_index)
    # defdomain the domain of definition of this potential (x_min,x_max,y_min,y_max),
    # num_states list of num of states in [x axis , y axis]

    x_min, x_max, y_min, y_max = defdomain

    # generate a grid of guesses dependent on the microstates to cover the whole potential surface
    guess = [
        np.linspace(x_min, x_max, int(num_states[0] / 2 + 1)),
        np.linspace(y_min, y_max, int(num_states[1] / 2 + 1)),
    ]
    guess = [[[i, j] for j in guess[1]] for i in guess[0]]
    min_list = []
    for i in guess:  # get min from each guess
        for j in i:
            min_list.append(fmin(potentialeval, j, disp=False))

    # sort out all the minima that may be calculated multiple times
    min_list = np.array(min_list).round(3)
    min_list = [list(i) for i in min_list]
    min_copy = min_list[::1]
    for i in sorted(list(range(len(min_list))), reverse=True):
        min_copy.pop(i)
        if min_list[i] in min_copy:
            min_list.pop(i)

    for i in range(len(min_list)):  # append the potential to each minimum
        popped = min_list.pop(0)
        popped.append(potentialeval(popped))
        min_list.append(popped)

    # define states
    states_h = [
        np.arange(x_min + 1 / (2 * num_states[0]) * (x_max - x_min), x_max, (x_max - x_min) / num_states[0]),
        np.arange(y_min + 1 / (2 * num_states[1]) * (y_max - y_min), y_max, (y_max - y_min) / num_states[1]),
    ]
    states_h = [[[j, i] for j in states_h[0]] for i in states_h[1]]
    states = []
    for i in states_h:
        for j in i:
            states.append(j)
    sts_xwidth = (x_max - y_min) / num_states[0] / 2
    sts_ywidth = (y_max - y_min) / num_states[1] / 2

    collect = []
    for st_idx in range(len(states)):

        # calculate the minimum in each state
        def help_fun(pos):
            if abs(pos[0] - states[st_idx][0]) < sts_xwidth and abs(pos[1] - states[st_idx][1]) < sts_ywidth:
                return potentialeval(pos)
            else:
                return math.inf

        st_min, st_minval, wfa, wfb, wfc = fmin(help_fun, states[st_idx], disp=False, full_output=True)

        # every state get"s assigned to its nearest local minimum
        dist = [(st_min[0] - lmin[0]) ** 2 + (st_min[1] - lmin[1]) ** 2 for lmin in min_list]
        for i in range(len(dist)):
            if dist[i] == min(dist):
                min_belong = i

        limit = min_list[min_belong][2] + threshhold

        if st_minval < limit:
            collect.append((min_belong, st_idx))

    ## Print output
    st_lst = []
    for i in range(len(min_list)):
        st = []
        for j in collect:
            if j[0] == i:
                st.append(j[1])
        if disp == True:
            print("The following states belong to the minimum at ({},{}):  ".format(min_list[i][0], min_list[i][1]), st)
        st_lst.append(st)
    return st_lst


##### calculate the first passage time distribution #####
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

    p = get_stat_distribution(T)
    fp_full = np.zeros(length)

    # weight by prob that TRAJECTORY starts in start
    weight = np.zeros(n)
    for j in start:
        for k in range(n):
            if k not in start:
                weight[j] += p[k] * T[k, j]

    weight = weight / sum(weight)  # normalise
    for j in range(N_start):
        fp_full += weight[start[j]] * fp[j]

    fp_full = fp_full / sum(fp_full)  # corrects for very small numercial error

    return fp_full


############################################

# ENTRPOY PRODUCTION

############################################


##### calculates the local entropy production from the model as well as the theoretical value from state 'start' to state 'end' #####
def entroprod(data_index, Tij, start, end):

    potentialeval, defdomain, num_states = Gauss_potential(data_index)  # get the potential function
    x_min, x_max, y_min, y_max = defdomain

    # defines microstates, 'states' describes a list of all microstates defined by their center
    states_h = [
        np.arange(x_min + 1 / (2 * num_states[0]) * (x_max - x_min), x_max, (x_max - x_min) / num_states[0]),
        np.arange(y_min + 1 / (2 * num_states[1]) * (y_max - y_min), y_max, (y_max - y_min) / num_states[1]),
    ]
    states_h = [[[j, i] for j in states_h[0]] for i in states_h[1]]
    states = []
    for i in states_h:
        for j in i:
            states.append(j)

    f = get_force(data_index)  # external force strength (only works for a constant force)
    T = 5  # need to be careful here if you change the temperature T

    # calculates the average potential within a micro state
    def avpot(pos):
        xwidth = x_max / num_states[0]
        ywidth = y_max / num_states[1]
        x_pos = np.linspace(-xwidth, xwidth, 10)
        x_pos += pos[0]
        y_pos = np.linspace(-ywidth, ywidth, 10)
        y_pos += pos[1]
        guess = np.array([[i, j] for j in y_pos for i in x_pos])
        pot_lst = []
        for i in guess:
            pot_lst.append(potentialeval(i))
        return np.mean(pot_lst)

    # calculate distance between the two states, accounts for the direction of the transition in question
    val = states[start][0] - states[end][0]
    if val > x_max / 2:
        dist = val - x_max
    elif val < -x_max / 2:
        dist = x_max + val
    elif val == x_max / 2 or val == -x_max / 2:
        dist = 0
    else:
        dist = val

    # calculate theoretical local entropy production as well as the one calculated from the model
    theo_val = (-f * dist + avpot(states[start]) - avpot(states[end])) / T
    cutoff = 0.0001  # avoid transitions which are not sufficiently sampled

    if Tij[start][end] and Tij[end][start] > cutoff:
        meas_val = np.log((Tij[start][end] / Tij[end][start]))
        return theo_val, meas_val

    else:
        return theo_val, None


##### caluculated global entropy production for a given transition matrix #####
def glob_entropy_prod(transition_matrix):
    Tij = transition_matrix
    pi = get_stat_distribution(Tij)
    entprod = []
    for i in range(len(Tij)):
        for j in range(len(Tij[i])):
            if Tij[i][j] and Tij[j][i] > 0:
                entprod.append(pi[i] * Tij[i][j] * np.log(pi[i] * Tij[i][j] / pi[j] / Tij[j][i]))
    return sum(entprod)


############################################

# REWEIGHTING ANALYSIS

############################################


##### plot the stationary distributions (for a reference, traget, and reweighted model) in a certain range #####
def plot_stat_distr(
    ref_data, target_data, lagtime, states
):  # states = [from,to] describes the range in terms of the microstates

    transition_matrix = get_transition_matrix(target_data, lagtime)
    vl1 = get_stat_distribution(transition_matrix)

    transition_matrix = get_transition_matrix(ref_data, lagtime)
    vl2 = get_stat_distribution(transition_matrix)

    transition_matrix = reweighting(ref_data, target_data, lagtime)
    vl3 = get_stat_distribution(transition_matrix)

    fig, ax = plt.subplots(dpi=150)
    ax.plot(np.arange(0, len(vl2[states[0] : states[1]]), 1), vl2[states[0] : states[1]], label="sim ref")
    ax.plot(np.arange(0, len(vl1[states[0] : states[1]]), 1), vl1[states[0] : states[1]], label=r"sim target")
    ax.plot(np.arange(0, len(vl3[states[0] : states[1]]), 1), vl3[states[0] : states[1]], label="rew")
    ax.legend(loc="best")
    ax.set_xlabel("No. microstate on x axis")
    plt.show()


##### plot the difference of the local entropy production between theory and model from a starting state #####
def plot_local_entprod(ref_data, target_data, lagtime, start):

    shape = get_shape(ref_data)

    end_lst = np.arange(0, shape[0] * shape[1])

    # calculate difference between theoretical value and reweigthed value for every transition from the starting state
    Tij = reweighting(ref_data, target_data, lagtime)
    tlst = []
    mlst = []
    for end in end_lst:
        theo_val, meas_val = entroprod(target_data, Tij, start, end)
        tlst.append(theo_val)
        mlst.append(meas_val)
    diff_rew = np.zeros(shape[0] * shape[1])
    for i in range(len(tlst)):
        if mlst[i] == None:
            diff_rew[i] = None
        else:
            diff_rew[i] = -tlst[i] + mlst[i]
    diff_rew = diff_rew.reshape(shape[1], shape[0])  # reshape for plotting

    # same for the difference of theretical value and simulation model value
    Tij = get_transition_matrix(target_data, lagtime)
    tlst = []
    mlst = []
    for end in end_lst:
        theo_val, meas_val = entroprod(target_data, Tij, start, end)
        tlst.append(theo_val)
        mlst.append(meas_val)
    diff_sim = np.zeros(shape[0] * shape[1])
    for i in range(len(tlst)):
        if mlst[i] == None:
            diff_sim[i] = None
        else:
            diff_sim[i] = -tlst[i] + mlst[i]
    diff_sim = diff_sim.reshape(shape[1], shape[0])

    # plotting
    fig, ax = plt.subplots(dpi=200)
    p = ax.imshow(diff_rew, interpolation="nearest", origin="lower", cmap="coolwarm", vmin=-0.5, vmax=0.5)

    fig.colorbar(p, shrink=0.3, extend="both")
    ax.plot(
        np.mod(start, shape[0]), int(start / shape[0]), marker="x", markersize=5, linestyle="", label="starting point"
    )
    ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("No. of micro-state on x-axis", fontsize=8)
    ax.set_ylabel("No. of micro-state on y-axis", fontsize=8)
    ax.tick_params(axis="both", labelsize=8)
    ax.set_title("reweighted model", fontsize=8)
    plt.show()

    fig, ax = plt.subplots(dpi=200)
    p = ax.imshow(diff_sim, interpolation="nearest", origin="lower", cmap="coolwarm", vmin=-0.5, vmax=0.5)

    ax.plot(
        np.mod(start, shape[0]), int(start / shape[0]), marker="x", markersize=5, linestyle="", label="starting point"
    )
    ax.legend(loc="best", fontsize=8)
    fig.colorbar(p, shrink=0.3, extend="both")
    ax.set_xlabel("No. of micro-state on x-axis", fontsize=8)
    ax.set_ylabel("No. of micro-state on y-axis", fontsize=8)
    ax.tick_params(axis="both", labelsize=8)
    ax.set_title("simulated model", fontsize=8)
    plt.show()


##### plot the mean first passage times for simulation and reweighted models #####
def plot_MFPTs(
    ref_data, target_lst, lagtime, length, start, end
):  # target_lst describes a list of all targets, start and end are the indices of the different wells

    min_states = state_assignement(Gauss_potential, ref_data, threshhold=0.01, disp=False)

    fig, ax = plt.subplots(dpi=150)
    x = np.arange(length)

    # calculate MFPTs from simulation
    mfptsim_lst = []
    for target_data in target_lst:
        ysim1 = first_passage(min_states[start], min_states[end], get_transition_matrix(target_data, lagtime), length)
        mfptsim_lst.append(np.dot(x, ysim1))
    ax.plot(np.arange(len(target_lst)), mfptsim_lst, linestyle="", marker="+", label="simulation")

    # same for reweighting
    mfptrew_lst = []
    for target_data in target_lst:
        yrew = first_passage(min_states[start], min_states[end], reweighting(ref_data, target_data, lagtime), length)
        mfptrew_lst.append(np.dot(x, yrew))
    ax.plot(np.arange(len(target_lst)), mfptrew_lst, linestyle="", marker="+", label="reweighting")

    # plotting
    ax.set_xlabel("extf", fontsize=15)
    ax.set_ylabel("MFPT", fontsize=15)
    ax.legend(loc="best")
    plt.show()


##### plot the mean first passage times for simulation and reweighted models with errorbars #####
def plot_MFPTs_error(
    ref_data, target_lst, lagtime, length, start, end, N, figax=False
):  # target_lst describes a list of all targets, start and end are the indices of the different wells, N is the amount of data used for the calculation

    min_states = state_assignement(Gauss_potential, ref_data, threshhold=0.01)
    start = min_states[start]
    end = min_states[end]

    x = np.arange(length)

    # calculate MFPTs for all simulation and reweighted models
    mfpt_sim_lst = []
    mfpt_rew_lst = []
    for target_data in target_lst:
        for i in range(0, N):
            target_data2 = str(target_data) + "_{}".format(i)
            ref_data2 = str(ref_data) + "_{}".format(i)

            Tij_sim = get_transition_matrix(target_data2, lagtime)
            Tij_rew = reweighting(ref_data2, target_data2, lagtime)
            mfpt_sim_lst.append(np.dot(x, first_passage(start, end, Tij_sim, length)))
            mfpt_rew_lst.append(np.dot(x, first_passage(start, end, Tij_rew, length)))

    # take the mean and error of these values
    mu_sim = []
    err_sim = []
    for i in range(len(target_lst)):
        mu_sim.append(np.mean(mfpt_sim_lst[i * N : N * i + N]))
        err_sim.append(np.std(mfpt_sim_lst[i * N : N * i + N]) / np.sqrt(N))

    mu_rew = []
    err_rew = []
    for i in range(len(target_lst)):
        mu_rew.append(np.mean(mfpt_rew_lst[i * N : N * i + N]))
        err_rew.append(np.std(mfpt_rew_lst[i * N : N * i + N]) / np.sqrt(N))

    # plotting
    ax.errorbar(
        np.arange(len(targetlst)),
        mu_sim,
        yerr=err_sim,
        linestyle="",
        marker=".",
        label="simulation",
        capsize=2,
        markersize=1,
    )
    ax.errorbar(
        np.arange(len(target_lst)),
        mu_rew,
        yerr=err_rew,
        linestyle="--",
        marker="x",
        label="reweighting",
        capsize=0,
        markersize=1,
    )

    ax.set_xlabel("$f_{ext}$")
    ax.set_ylabel("MFPT [\tau]")
    ax.legend(loc="best")
    plt.show()
    if figax == True:
        return fig, ax


##### plot the global entropy production for simulation and reweighted models #####
def plot_glob_entprod(ref_data, target_lst, lagtime):  # target_lst describes a list of all targets

    fig, ax = plt.subplots(dpi=150)

    # calculate glob_entprod from simulation
    entsim_lst = []
    for target_data in target_lst:
        entsim_lst.append(glob_entropy_prod(get_transition_matrix(target_data, lagtime)))

    ax.plot(np.arange(len(target_lst)), entsim_lst, linestyle="", marker="+", label="simulation")

    # same for reweighting
    entrew_lst = []
    for target_data in target_lst:
        entrew_lst.append(glob_entropy_prod(reweighting(ref_data, target_data, lagtime)))
    ax.plot(np.arange(len(target_lst)), entrew_lst, linestyle="", marker="+", label="reweighting")

    # plotting
    ax.set_xlabel("extf", fontsize=15)
    ax.set_ylabel("glob entropy prod", fontsize=15)
    ax.legend(loc="best")
    plt.show()


##### plot the global entropy production for simulation and reweighted models with errorbars#####
def plot_globalent_error(
    ref_data, target_data_lst, lagtime, N, figax=False
):  # N is the amount of data used for the calculation

    # calculate simulation and reweighted values for all data
    ent_sim_lst = []
    ent_rew_lst = []
    for target_data in target_data_lst:
        for i in range(0, N):
            target_data2 = str(target_data) + "_{}".format(i)
            ref_data2 = str(ref_data) + "_{}".format(i)
            Tij_sim = get_transition_matrix(target_data2, lagtime)
            Tij_rew = reweighting(ref_data2, target_data2, lagtime)
            ent_sim_lst.append(glob_entropy_prod(Tij_sim))
            ent_rew_lst.append(glob_entropy_prod(Tij_rew))

    # take the mean and error of these values
    mu_sim = []
    err_sim = []
    for i in range(len(target_data_lst)):
        mu_sim.append(np.mean(ent_sim_lst[i * N : N * i + N]))
        err_sim.append(np.std(ent_sim_lst[i * N : N * i + N]) / np.sqrt(N))
    mu_rew = []
    err_rew = []
    for i in range(len(target_data_lst)):
        mu_rew.append(np.mean(ent_rew_lst[i * N : N * i + N]))
        err_rew.append(np.std(ent_rew_lst[i * N : N * i + N]) / np.sqrt(N))

    # plotting
    fig, ax = plt.subplots(dpi=150)
    ax.plot(np.arange(len(target_data_lst)), mu_sim, linestyle="", marker=".", label="simulation")
    ax.plot(np.arange(len(target_data_lst)), mu_rew, linestyle="--", label="reweighting")
    ax.plot(np.arange(len(target_data_lst))[0], mu_rew[0], linestyle="", marker=".", markersize=15)
    ax.set_xlabel("extf")
    ax.set_ylabel("Entrpoy production")
    ax.legend(loc="best")
    plt.show()

    if figax == True:
        return fig, ax


##### plot FPTD of the simulated and reweighted model #####
def plot_FPTD(
    ref_data, target_data, lagtime, length, start, end, figax=False
):  # start and end are the indices of the different wells

    min_states = state_assignement(Gauss_potential, ref_data, threshhold=0.01, disp=False)

    fig, ax = plt.subplots(dpi=150)
    x = np.arange(length)
    y = first_passage(min_states[start], min_states[end], get_transition_matrix(target_data, lagtime), length)
    ax.plot(x, y, label="simulation (extf = 9)", marker="x", markersize=2, linestyle="")

    y = first_passage(min_states[start], min_states[end], reweighting(ref_data, target_data, lagtime), length)
    ax.plot(x, y, label="reweighting from extf = {}".format(str(ref_data)[1]), marker="x", markersize=2, linestyle="")

    ax.set_xlabel("#steps", fontsize=15)
    ax.set_ylabel("p", fontsize=15)
    ax.legend(loc="best")
    plt.show()
    if figax == True:
        return fig, ax


##### caluculated FPT moments beforehand and save them in a file #####
def calc_fpt_moments(calc_lst, lagtime, length=1000):  # calc_lst = data_set in the form of ['7x06']

    x = np.arange(length)
    idx = np.arange(0, 10, 1)  # needs to be adjusted for sets with a different form

    for data in calc_lst:
        data_lst = [int(data.replace("x", "%s" % i)) for i in idx]
        min_states = state_assignement(Gauss_potential, data_lst[0], threshhold=0.01, disp=False)

        # account for all different transitions between wells
        for ref_data in data_lst:
            for start_idx in range(len(min_states)):
                for end_idx in range(len(min_states)):
                    if start_idx != end_idx:

                        # data path needs to be adjusted
                        with open(
                            "../data/FPT/FPT_{}_{}{}.txt".format(ref_data, start_idx, end_idx), mode="w+"
                        ) as file:

                            start = min_states[start_idx]
                            end = min_states[end_idx]

                            for target_data in data_lst:
                                Tij_sim = get_transition_matrix(target_data, lagtime)
                                Tij_rew = reweighting(ref_data, target_data, lagtime)
                                fpt_sim = first_passage(start, end, Tij_sim, length)
                                fpt_rew = first_passage(start, end, Tij_rew, length)

                                # calculate first three moments
                                mu_fpt_sim = np.dot(x, fpt_sim)
                                mu_fpt_rew = np.dot(x, fpt_rew)
                                sig_fpt_sim = np.sqrt(np.dot(x**2, fpt_sim) - mu_fpt_sim**2)
                                sig_fpt_rew = np.sqrt(np.dot(x**2, fpt_rew) - mu_fpt_rew**2)
                                skew_fpt_sim = (
                                    np.dot(x**3, fpt_sim) - 3 * sig_fpt_sim**2 * mu_fpt_sim - mu_fpt_sim**3
                                ) / sig_fpt_sim**3
                                skew_fpt_rew = (
                                    np.dot(x**3, fpt_rew) - 3 * sig_fpt_rew**2 * mu_fpt_rew - mu_fpt_rew**3
                                ) / sig_fpt_rew**3

                                # save in file
                                file.write(
                                    "{}, {}, {}, {}, {}, {} \n".format(
                                        mu_fpt_sim, mu_fpt_rew, sig_fpt_sim, sig_fpt_rew, skew_fpt_sim, skew_fpt_rew
                                    )
                                )


##### plot the FPT moments if calculated beforehand #####
def plot_FPT_moments(data_index, start, end, MFPT=True, sig=False, skew=False):
    file = "/home/theorie/wittenstein/data/FPT/FPT_{}_{}{}.txt".format(
        data_index, start, end
    )  # data path needs to be adjusted
    data = np.loadtxt(file, delimiter=",")

    if MFPT == True:
        fig, ax = plt.subplots(dpi=150)
        ax.plot(np.arange(len(data)), data[:, 0], linestyle="", marker="+", label="simulation")
        ax.plot(np.arange(len(data)), data[:, 1], linestyle="", marker="+", label="reweighting")
        ax.set_xlabel("extf", fontsize=15)
        ax.set_ylabel("MFPT", fontsize=15)
        ax.legend(loc="best")
        plt.show()

    if sig == True:
        fig, ax = plt.subplots(dpi=150)
        ax.plot(np.arange(len(data)), data[:, 2], linestyle="", marker="+", label="simulation")
        ax.plot(np.arange(len(data)), data[:, 3], linestyle="", marker="+", label="reweighting")
        ax.set_xlabel("extf", fontsize=15)
        ax.set_ylabel("sig", fontsize=15)
        ax.legend(loc="best")
        plt.show()

    if skew == True:
        fig, ax = plt.subplots(dpi=150)
        ax.plot(np.arange(len(data)), data[:, 4], linestyle="", marker="+", label="simulation")
        ax.plot(np.arange(len(data)), data[:, 5], linestyle="", marker="+", label="reweighting")
        ax.set_xlabel("extf", fontsize=15)
        ax.set_ylabel("skewness", fontsize=15)
        ax.legend(loc="best")
        plt.show()
