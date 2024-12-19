import numpy as np


def define_cross(minpos,cut,dist,ms): 
	states = len(minpos)
	cross = np.zeros((states,2))
	if (cut[0] <  minpos[0]):
		for k in range(2,states):
			cross[k,0] = min(minpos[k-2]+dist-1,minpos[k]-dist+1)
			cross[k,1] = max(minpos[k-2]+dist,minpos[k]-dist)
		
		cross[0,1] = min(minpos[0]-dist,minpos[1]+dist)
		cross[0,0] = max(minpos[0]-dist+1,minpos[1]+dist-1)
		cross[1,1] = min(minpos[1]-dist,minpos[states-1]+dist)
		cross[1,0] = max(minpos[1]-dist+1,minpos[states-1]+dist-1)
		
		outpos = np.where(cross >= ms)
		cross[outpos] -= ms

		outpos = np.where(cross < 0)
		cross[outpos] += ms

	else:	
		print("TODO : case cut[0] > minpos[0] ")			
	
#	print(cross)
	return cross	


