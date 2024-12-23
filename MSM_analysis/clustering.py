import numpy as np
import os
from sklearn.cluster import KMeans


def dist(p1, p2):
	box = [3.,1.]  #adjust boxsize 
	x = np.minimum( np.abs(p1[0] - p2[0]), np.minimum(p1[0],p2[0]) +box[0] - np.maximum(p1[0],p2[0])  )
	y = np.minimum( np.abs(p1[1] - p2[1]), np.minimum(p1[1],p2[1]) +box[1] - np.maximum(p1[1],p2[1])  )
	
	return np.linalg.norm(x)


dimstart = 1
dimstop = 3
dim = dimstop - dimstart
lag = 10
cl = 300

path = "/data/pckr194/bause/single_particle/trajectory_8000_0.dat"
tra_len = int(os.popen("wc -l "+path).readline().split()[0])
fi =  open(path)
N = int(tra_len / lag )

print(N)
#data = np.loadtxt("/data/pckr194/bause/single_particle/trajectory_8000_0.dat")
#d = data[::lag,dimstart:dimstop]
d = np.zeros((N,2))
for i in range(0,N):
	line= fi.readline()	
	if(i % lag == 0):
		if line.strip():
			splt = line.split()
			d[i,0] = float(splt[1])
			d[i,1] = float(splt[2])



KMeans.euclidean_distances = dist

clusters = KMeans(n_clusters = cl,precompute_distances=False,algorithm='full', verbose =1, n_init=5,tol = 1e-04).fit(d)



np.savetxt("test.dat",clusters.cluster_centers_)


np.savetxt("/data/isilon/bause/single_particle/cluster_8000.dat",clusters.cluster_centers_)



