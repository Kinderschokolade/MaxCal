import numpy as np

def autocorr(x,maxpos):
	result = np.correlate(x, x, mode = 'full')
	#maxcorr = np.argmax(result)
	#print 'maximum = ', result[maxcorr]
	#result = result / result[maxcorr]
	return result[len(x)-1:len(x)-1+maxpos]


def load_column(filename, column):
	p = []
	i = 0
	with open(filename) as inf:
		for line in inf:
			parts = line.split() # split line into parts
			p.append(float(parts[column]))
	p=np.array(p)
	return(p)

p = []
corr = []
column = 4
maxpos = 100
name = "/data/isilon/bause/single_particle/trajectory_1405.dat"
p = load_column(name,column)
p = p[0:100000]
print(len(p)) 
pav = np.average(p)
print(pav)
corr = autocorr(p-pav,maxpos)
f = 'corr_1405.dat'
np.savetxt(f, corr, delimiter='\t')
D = np.trapz(corr,dx= 0.01)
print(D)
print(D/pav) #actually *f , whatever that is  








