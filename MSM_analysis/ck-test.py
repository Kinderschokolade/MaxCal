import numpy as np

def escape_tr(L,q,c,relax_pos,lvl,ms):
	pi = np.zeros((lvl,ms))
	for j in range(lvl):
		for i in relax_pos[j]:
			pi[j,i] = q[i]
		pi[j] = pi[j] / sum(pi[j,:]) # starting vector

	p = np.zeros(2*lvl)
	px = np.zeros((ms,lvl))
	
	for k in range(lvl):
		# transition matrix from  ms to A
		for l in relax_pos[k]:
			for m in relax_pos[k]:
				px[l,k] += L[l,m]
		# transition matrix from A to A
		for l in relax_pos[k]:	
			p[k] += pi[k,l]* px[l,k] 
	
		# err estmation from A to A
		p[lvl+k] = np.sqrt(k*(p[k]-p[k]*p[k]) / sum(sum(c[y*ms+z] for z in range(ms)) for y in relax_pos[k]))

	return p
		



def escape(pi,T,relax_pos,lvl,steps,lag):	
	psum = np.zeros((lvl+1,steps))
	ms = pi.size
	p = np.zeros((lvl,ms))
	p_old = np.zeros((lvl,ms))
	for k in range(lvl):
		for i in relax_pos[k]:
			psum[k+1,0] += pi[i]
	
	for j in range(lvl):
		for i in relax_pos[j]:
			p_old[j,i] = pi[i]
		p_old[j] = p_old [j] / sum(p_old[j,:]) # starting vector
		psum[j+1,0] = 1.

	psum[0,0] = 0


	for s in range(1,steps):
		# markov model estimation
		for k in range(lvl):
			for i in range(ms):
				for j in range(ms):
					p[k,j] += p_old[k,i] * T[i,j] 

		for k in range(lvl):
			for i in relax_pos[k]:
				psum[k+1,s] += p[k,i]  # prob to stay in same basin
		psum[0,s] = s*lag	
		p_old[:,:] = p[:,:]
		p[:,:] =0.
	

	return psum


ident = '8000'
path0  = "/data/isilon/bause/"

cols = 999

data = np.loadtxt(path0+"single_particle/count_"+ str(ident) +".dat", dtype='i',usecols = (range(1,cols+1)))

ms = 300
tau = 800

lvl =3

c = data[0:ms*ms]
p = data[ms*ms:ms*ms+ms]   

L = np.zeros((ms,ms))


steps = 100;
p_ck = np.zeros((lvl+1,steps))
p_tr = np.zeros((2*lvl+1,c.shape[1]))


relax_pos = np.zeros((lvl,4))
relax_pos[0] = [124, 125, 154, 155]
relax_pos[1] = [134, 135, 164, 165]
relax_pos[2] = [144, 145, 174, 175]


for i in range(c.shape[1]):
	T = np.transpose(c[:,i].reshape((ms,ms)))
	q = p[:,i]
	
	for k in range(ms):
		for l in range(ms):
			if (q[k] > 0.):
				L[k,l] = T[k,l] / q[k]	
			else:
				L[k,l] = 0.
				L[k,k] = 1. 

	p_tr[0,i] = i+1;
	p_tr[1:,i] = escape_tr(L,q,c[:,i],relax_pos,lvl,ms)
	if (i == tau):	
		L_check = L
		#L = construct_MSM(L,F,Sp,ms)
		T_check = T

		np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/ps_'+str(ident)+'.dat',q/sum(q))
		p_check = q/sum(q)
		#escapetime for MSM chapman kolmogorov test
		p_ck = escape(p_check,L,relax_pos,lvl,steps,tau)

		np.savetxt(path0+'single_particle/'+str(int(tau))+'/MFPT/cktest_'+str(ident)+'.dat',np.transpose(p_ck))



np.savetxt(path0+"single_particle/lagtime/ck_"+ident+".dat",np.transpose(p_tr))	







