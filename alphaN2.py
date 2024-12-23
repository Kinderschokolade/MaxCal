import argparse
import subprocess
import numpy as np
import multiprocessing as mp

parser = argparse.ArgumentParser()
arg_names = ['output', 'T1', 'T2', 'T3', 'U1', 'U2', 'U3', 'gamma', 'dT', 'fullT', 'k']
for arg in arg_names:
    parser.add_argument(f'-{arg[0]}', f'--{arg}', help=f'{arg} value')

args = parser.parse_args()
path = '/data/isilon/bause/single_particle/'

def load_column(filename, column):
    with open(filename) as inf:
        return np.array([float(line.split()[column]) for line in inf])

def load_row(filename, grep):
    thing = subprocess.Popen(f"grep {grep} {filename}", shell=True, stdout=subprocess.PIPE).communicate()[0]
    out = thing.split()[1:]
    return np.asarray(out, dtype=np.float32)

def run_Sim(U1, U2, U3, T1, T2, T3, gamma, argsout, fullT, potk, dT, filepos, i):
    U1vary = str(int(U1) + i * 0.5 + 0.5)
    filename = f"{int(argsout) + filepos}{'0' if i < 0 else ''}{i + 10}"
    cmd = f"~/code/single_particle/a.out -U1 {U1vary} -U2 {U2} -U3 {U3} -T1 {T1} -T2 {T2} -T3 {T3} -gamma {gamma} -o {filename} -fullT {fullT} -pot 2 -k {potk} -dT {dT} -it {i + 10}"
    print(cmd)
    sim = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    sim.wait()

def run_sim_u1(uu11):
    return run_Sim(uu11, args.U2, args.U3, args.T1, args.T2, args.T3, args.gamma, args.output, args.fullT, args.k, args.dT, 123, 0)


# run 
if __name__ == '__main__':
    pool = mp.Pool()
    run_sim_u1_range = range(-int(args.U1), int(args.U2))
    pool.map_async(run_sim_u1, run_sim_u1_range).get(9999999)
