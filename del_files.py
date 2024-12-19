import shutil 
 
 
for i in range(4,300): 
        if (i%10 != 0): 
                path = "/data/isilon/bause/single_particle/"+str(i) 
                #print(path) 
                shutil.rmtree(path)

