import analyse_2D as ana 

ident = 9000
odent = list(range(9100,9800,200))
#odent = 8900

lagtime = 200
length = 1000


states= ana.state_assignement(ana.Gauss_potential,ident)

start = 1
end = 2

#ana.plot_FPTD(ident,odent,lagtime,length,start,end,figax=True)
ana.plot_MFPTs(ident,odent,lagtime,length,start,end)



