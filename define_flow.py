import numpy as np
	
def pos_pot(i,lvl,cut):
	pos = lvl
	for k in range(lvl):
		if (i <= cut[k]):
			pos = k
			break
	return pos


def def_flow_direction(ms):
	r = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
			r[i,j] = i -j

	return r


def define_flow_direction(ms,cut): # works much easier using a dictionary... whatever
	r = np.zeros((ms,ms))
	lvl = np.size(cut)
	
	for i in range(ms):
		pos = pos_pot(i,lvl,cut)
		for j in range(ms):
			posto = pos_pot(j,lvl,cut)	
			if (pos == posto):
#				r[i,j] = np.sign(j-i) 
				r[i,j] = j-i
			if (pos + 1 ==  posto):
#				r[i,j] = 1
				r[i,j] = j-i
			if (pos  == posto +1 ):
#				r[i,j] = -1
				r[i,j] = j-i

			if (posto == 0 and pos == lvl):
#				r[i,j] = 1
				r[i,j] = ms + j - i
			if (posto == lvl and pos == 0):
#				r[i,j] = -1
				r[i,j] = -i -ms +j 
			if (posto == lvl-1 and pos == 0):
#				r[i,j] = -1
				r[i,j] = -i -ms +j
			if (posto == 0 and pos == lvl-1): #lvl-1 to lvl is covered
#				r[i,j] = 1
				r[i,j] = ms +j -i	
			if (pos == 1 and posto == lvl):
#				r[i,j] = -1
				r[i,j] = -i-ms+j
			if (pos == lvl and posto == 1):
#				r[i,j] = 1
				r[i,j] = ms+j-i
	return r

def read_cut(ident):
	xc = []
	U = []
	for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
		cont = line.split()
		if(cont[0][:2] == 'xc'):
			xc.append (float(cont[1]))
		if (cont[0] == 'microstates'):
			ms = int(cont[1])
		if (cont[0][0] == 'U'):
			U.append(float(cont[1]))


	reg = []
	N = len(xc)
	for i in range(1, N):
		reg.append(int(round(xc[i] *ms)))

	if (reg[N-2] >= ms):
		reg.insert(0, int(round(ms * xc[0])))
		reg.pop(N-1)
	cut = []
	if (U[0] > U[1]):
		cut.append(round((ms - reg[N-2] + reg[0]) /2))
		for i in range(2,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ))
	else:
		for i in range(0,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ))
	
	return ms,cut

def read_cut(path, ident):
	xc = []
	U = []
	for line in open(path+"single_particle/param_"+str(ident)+".dat","r"):
		cont = line.split()
		if(cont[0][:2] == 'xc'):
			xc.append (float(cont[1]))
		if (cont[0] == 'microstates'):
			ms = int(cont[1])
		if (cont[0][0] == 'U'):
			U.append(float(cont[1]))


	reg = []
	N = len(xc)
	for i in range(1, N):
		reg.append(int(round(xc[i] *ms)))

	if (reg[N-2] >= ms):
		reg.insert(0, int(round(ms * xc[0])))
		reg.pop(N-1)
	cut = []
	if (U[0] > U[1]):
		cut.append(round((ms - reg[N-2] + reg[0]) /2))
		for i in range(2,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ))
	else:
		for i in range(0,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ))
	
	return ms,cut



def read_minpos(ident):
	xc = []
	U = []
	for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
		cont = line.split()
		if(cont[0][:2] == 'xc'):
			xc.append (float(cont[1]))
		if (cont[0] == 'microstates'):
			ms = int(cont[1])
		if (cont[0][0] == 'U'):
			U.append(float(cont[1]))


	reg = []
	N = len(xc)
	for i in range(1, N):
		reg.append(int(round(xc[i] *ms)))

	if (reg[N-2] > ms):
		reg.insert(0, int(round(ms * xc[0])))
		reg.pop(N-1)


	cut = []
	if (U[0] < U[1]):
		cut.append(round((ms - reg[N-2] + reg[0]) /2) -1)
		for i in range(1,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ) -1)
	else:
		for i in range(0,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ) -1)

	return ms,cut



def read_minpos(path, ident):
	xc = []
	U = []
	for line in open(path+"single_particle/param_"+str(ident)+".dat","r"):
		cont = line.split()
		if(cont[0][:2] == 'xc'):
			xc.append (float(cont[1]))
		if (cont[0] == 'microstates'):
			ms = int(cont[1])
		if (cont[0][0] == 'U'):
			U.append(float(cont[1]))


	reg = []
	N = len(xc)
	for i in range(1, N):
		reg.append(int(round(xc[i] *ms)))

	if (reg[N-2] > ms):
		reg.insert(0, int(round(ms * xc[0])))
		reg.pop(N-1)


	cut = []
	if (U[0] < U[1]):
		cut.append(round((ms - reg[N-2] + reg[0]) /2) -1)
		for i in range(1,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ) -1)
	else:
		for i in range(0,N-2,2):
			cut.append(round( (reg[i] + reg[i+1])/2 ) -1)

	return ms,cut

#
#def read_minpos(potential,box,Nx): # include potential fct, nD
#	dim= box.size()
#	x0 = np.zeros(dim)
#	bnds = ((0,box[0]),(0,box[1]))
#	minlist= [];
#	flag = 1;
#	tol = 1e-05
#	for i in range(Nx):
#		x0[0] = i/box[0]
##		for j in range(Ny):
##			x0[1] = j/box[1]
#			res = minimize(potential,x0,args=(),bounds=bnds, method='L-BFGS-B')
#			for k in range(len(minlist)):
#					if (tol > abs(res-minlist[k])):
#						flag = 0
#			
#			if flag:
#				minlist.append(res)
#			flag = 1
#
#	print(minlist)
#	return minlist
#









