from define_cross import define_cross
from define_flow import pos_pot
import numpy as np
import math

def get_pos(x,ms):
	pos = math.floor(x*ms)
	if(pos == ms):
		 pos -= 1
	return pos


def mfpt_trajectory_cts_reg(Gamma, blocks, minpos, cut, lvl, rancut,ms):
	pos = Gamma[0]
	pos_old = pos
	cross = define_cross(minpos,cut,rancut,ms)
	datalen_block = int(len(Gamma) / blocks)
	reg1 = pos_pot(pos_old,lvl,cut)
	if (reg1 == lvl):
		reg =0
	else:
		reg = reg1

	steps = np.zeros((lvl,lvl))
	regmax = lvl -1
	MFPT = np.zeros((blocks,lvl,lvl))
	MFPT_steps = np.zeros((blocks,lvl,lvl))
	for b in range(0,blocks):
		for j in range(0,datalen_block):
			if (abs(pos_old - pos) > ms/2):
				if(pos < pos_old):
					ran1 = range(-2,pos+1)
					ran2 = range(pos_old,ms+2)
				else:
					ran1 = range(-2,pos_old+1)
					ran2 = range(pos,ms+2)
			else:
				if (pos > pos_old):
					ran1 = range(pos_old,pos+1)
					ran2 = ran1
				else:
					ran1 = range(pos,pos_old+1)
					ran2 = ran1


			if (cross[reg,1] in ran1 or cross[reg,1] in ran2 or cross[reg,0] in ran1 or cross[reg,0] in ran2 ):
			# if new region was reached and minimum position was crossed	
				reg_old = reg
				if(cross[reg,0] in ran1 or cross[reg,0] in ran2):
					reg =  reg - 1 if reg > 0 else regmax
				else:
					reg = reg +1 if reg < regmax else 0

				#print(b*datalen_block+j,reg_old, reg, pos_old, pos, steps)

				MFPT[b,reg_old,reg] = (MFPT[b,reg_old,reg] * MFPT_steps[b,reg_old,reg] + steps[reg_old,reg]) / (MFPT_steps[b,reg_old,reg] +1.)
				MFPT_steps[b,reg_old,reg] +=1
				steps[:,:] +=1
				#steps[:,:] = 1
				steps[reg_old,reg] = 1
			else:
				steps[:,:] += 1

			pos_old = pos
			pos = Gamma[b*datalen_block+j]

		#print(MFPT[b])
		
	
	MFPT_av = np.zeros((lvl,lvl))
	MFPT_err = np.zeros((lvl,lvl))
	MFPT_stepsav = np.zeros((lvl,lvl))

	for b in range(blocks):
		MFPT_av += MFPT[b,:,:]
		MFPT_stepsav += MFPT_steps[b,:,:]

	MFPT_av /= blocks

	for b in range(blocks):
		temp = MFPT_av - MFPT[b,:,:]
		MFPT_err += temp *temp
	MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )
	return MFPT_av, MFPT_err, MFPT_stepsav

def get_reg(pos, minpos, rancut):
	for i in range(len(minpos)):
		if (pos >= (minpos[i] - rancut) and pos <= (minpos[i] + rancut)):
			return i
	
	return None

def mfpt_trajectory_area(Gamma, blocks, minpos, cut, lvl, rancut, ms):
	pos = Gamma[0]
	pos_old = pos
	datalen_block = int(len(Gamma) / blocks)
	num_reg = len(minpos)

	steps = np.zeros((num_reg,num_reg))
	MFPT = np.zeros((blocks,num_reg,num_reg))
	MFPT_steps = np.zeros((blocks,num_reg,num_reg))
	
	maxstep = 1000
	MFPT_hist = np.zeros((num_reg*num_reg,maxstep))

	reg = get_reg(pos, minpos,rancut)
	checkpos = np.zeros((num_reg,num_reg))
	if reg is not None:
		checkpos[reg,:] = 1

	for b in range(0,blocks):
		for j in range(0,datalen_block):
			if reg is not None:
				checkpos[reg,:] = 1 

				for l in range(num_reg):
					if (checkpos[l,reg] == 1 and l != reg):
						MFPT[b,l,reg] = (MFPT[b,l,reg] * MFPT_steps[b,l,reg] + steps[l,reg]) / (MFPT_steps[b,l,reg] +1.)
						MFPT_steps[b,l,reg] +=1
						if (steps[l,reg] < maxstep):
							MFPT_hist[l*3+reg,steps[l,reg]] += 1
						steps[l,reg] = 0
						checkpos[l,reg] = 0
			
			for i in range(num_reg):
				for k in range(num_reg):
					if (checkpos[i,k] == 1):
						steps[i,k] +=1 

			pos = Gamma[b*datalen_block+j]
			reg = get_reg(pos,minpos,rancut)

	MFPT_av = np.zeros((num_reg,num_reg))
	MFPT_err = np.zeros((num_reg,num_reg))
	MFPT_stepsav = np.zeros((num_reg,num_reg))

	for b in range(blocks):
		MFPT_av += MFPT[b,:,:]
		MFPT_stepsav += MFPT_steps[b,:,:]

	MFPT_av /= blocks

	for b in range(blocks):
		temp = MFPT_av - MFPT[b,:,:]
		MFPT_err += temp *temp
	MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )
	return MFPT_av, MFPT_err, MFPT_stepsav, MFPT_hist



def mfpt_trajectory_ms(Gamma, blocks, minpos, cut, lvl, rancut,ms):
	pos = Gamma[0]
	pos_old = pos
	datalen_block = int(len(Gamma) / blocks)

	steps = np.zeros((ms,ms))
	MFPT = np.zeros((blocks,ms,ms))
	MFPT_steps = np.zeros((blocks,ms,ms))
	

	checkpos = np.zeros((ms,ms))
	checkpos[pos,:] = 1

	for b in range(0,blocks):
		for j in range(0,datalen_block):
			checkpos[pos,:] = 1 
			for l in range(ms):
				if (checkpos[l,pos] == 1 and l != pos):
					MFPT[b,l,pos] = (MFPT[b,l,pos] * MFPT_steps[b,l,pos] + steps[l,pos]) / (MFPT_steps[b,l,pos] +1.)
					MFPT_steps[b,l,pos] +=1
					steps[l,pos] = 0
					checkpos[l,pos] = 0
			
			for i in range(ms):
				for k in range(ms):
					if (checkpos[i,k] == 1):
						steps[i,k] +=1 

			pos = Gamma[b*datalen_block+j]

	MFPT_av = np.zeros((ms,ms))
	MFPT_err = np.zeros((ms,ms))
	MFPT_stepsav = np.zeros((ms,ms))

	for b in range(blocks):
		MFPT_av += MFPT[b,:,:]
		MFPT_stepsav += MFPT_steps[b,:,:]

	MFPT_av /= blocks

	for b in range(blocks):
		temp = MFPT_av - MFPT[b,:,:]
		MFPT_err += temp *temp
	MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )
	return MFPT_av, MFPT_err, MFPT_stepsav


def mfpt_trajectory_ms_cross(Gamma, blocks, minpos, cut, lvl, rancut,ms, r):
	pos = Gamma[0]
	pos_old = pos
	datalen_block = int(len(Gamma) / blocks)

	steps = np.zeros((ms,ms))
	MFPT = np.zeros((blocks,ms,ms))

	checkpos = np.zeros((ms,ms))
	checkpos[pos,:] =1

	MFPT_steps = np.zeros((blocks,ms,ms))
	for b in range(0,blocks):
		for j in range(0,datalen_block):
			#if (pos != pos_old):
			# if new region was reached and minimum position was crossed	
			if (r[pos_old,pos] == -1):
				if (pos > pos_old):
					poslist = list(range(pos_old+1, pos+1))
				else:
					poslist = list(range(pos_old+1,ms)) + list(range(pos+1))
				
				
				#if pos_old in poslist:
				#	poslist.remove(pos_old)
				checkpos[pos,:] = 1

				for l in range(ms):
					for k in poslist:
						if (checkpos[l,k] == 1 and l != k):
						#MFPT[b,pos_old,k] = (MFPT[b,pos_old,k] * MFPT_steps[b,pos_old,k] + steps[pos_old,k]) / (MFPT_steps[b,pos_old,k] +1.)
						#MFPT_steps[b,pos_old,k] +=1
						#steps[pos_old,k] = 0	
							MFPT[b,l,k] = (MFPT[b,l,k] * MFPT_steps[b,l,k] + steps[l,k]) / (MFPT_steps[b,l,k] +1.)
							MFPT_steps[b,l,k] +=1
							steps[l,k] = 0
							checkpos[l,k] = 0
			elif (r[pos_old,pos] == 1):
				if (pos_old > pos):
					poslist = list(range(pos, pos_old))
				else:
					poslist = list(range(pos,ms)) + list(range(pos_old))
			
				checkpos[pos,:]  =1
				for l in range(ms):	
					for k in poslist:
						if (checkpos[l,k] ==1 and l != k):
					#MFPT[b,pos_old,k] = (MFPT[b,pos_old,k] * MFPT_steps[b,pos_old,k] + steps[pos_old,k]) / (MFPT_steps[b,pos_old,k] +1.)
					#MFPT_steps[b,pos_old,k] +=1
					#steps[pos_old,k] = 0
							MFPT[b,l,k] = (MFPT[b,l,k] * MFPT_steps[b,l,k] + steps[l,k]) / (MFPT_steps[b,l,k] +1.)
							MFPT_steps[b,l,k] +=1
							steps[l,k] = 0
							checkpos[l,k] = 0
			for i in range(ms):
				for k in range(ms):
					if (checkpos[i,k] ==1):
						steps[i,k] += 1
			pos_old = pos
			pos = Gamma[b*datalen_block+j]

	MFPT_av = np.zeros((ms,ms))
	MFPT_err = np.zeros((ms,ms))
	MFPT_stepsav = np.zeros((ms,ms))

	for b in range(blocks):
		MFPT_av += MFPT[b,:,:]
		MFPT_stepsav += MFPT_steps[b,:,:]

	MFPT_av /= blocks

	for b in range(blocks):
		temp = MFPT_av - MFPT[b,:,:]
		MFPT_err += temp *temp
	MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )
	return MFPT_av, MFPT_err, MFPT_stepsav


def mfpt_trajectory_cts_cross(Gamma, blocks, minpos, cut, lvl, rancut,ms, r):
	pos = get_pos(Gamma[0],ms)
	pos_old = pos
	datalen_block = int(len(Gamma) / blocks)

	steps = np.zeros((ms,ms))
	MFPT = np.zeros((blocks,ms,ms))

	MFPT_steps = np.zeros((blocks,ms,ms))
	for b in range(0,blocks):
		for j in range(0,datalen_block):
			#if (pos != pos_old):
			# if new region was reached and minimum position was crossed	
			if (r[pos_old,pos] == -1):
				if (pos > pos_old):
					poslist = list(range(pos_old+1, pos+1))
				else:
					poslist = list(range(pos_old+1,ms)) + list(range(pos+1))
				
				if pos_old in poslist:
					poslist.remove(pos_old)
				for k in poslist:
					MFPT[b,pos_old,k] = (MFPT[b,pos_old,k] * MFPT_steps[b,pos_old,k] + steps[pos_old,k]) / (MFPT_steps[b,pos_old,k] +1.)
					MFPT_steps[b,pos_old,k] +=1
					steps[pos_old,k] = 0	

				steps[:,:] +=1
			
			elif (r[pos_old,pos] == 1):
				if (pos_old > pos):
					poslist = list(range(pos, pos_old))
				else:
					poslist = list(range(pos,ms)) + list(range(pos_old))
					
				if pos_old in poslist:
					print("! "  )
					poslist.remove(pos_old)

				for k in poslist:
					MFPT[b,pos_old,k] = (MFPT[b,pos_old,k] * MFPT_steps[b,pos_old,k] + steps[pos_old,k]) / (MFPT_steps[b,pos_old,k] +1.)
					MFPT_steps[b,pos_old,k] +=1
					steps[pos_old,k] = 0
				steps[:,:] +=1
			else:
				steps[:,:] += 1

			pos_old = pos
			pos = get_pos(Gamma[b*datalen_block+j],ms)

	MFPT_av = np.zeros((ms,ms))
	MFPT_err = np.zeros((ms,ms))
	MFPT_stepsav = np.zeros((ms,ms))

	for b in range(blocks):
		MFPT_av += MFPT[b,:,:]
		MFPT_stepsav += MFPT_steps[b,:,:]

	MFPT_av /= blocks

	for b in range(blocks):
		temp = MFPT_av - MFPT[b,:,:]
		MFPT_err += temp *temp
	MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )
	return MFPT_av, MFPT_err, MFPT_stepsav





