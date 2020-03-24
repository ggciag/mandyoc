import numpy as np
import matplotlib.pyplot as plt

with open("param_1.5.3_2D.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx,Nz = int(line[0]),int(line[1])
	line = f.readline()
	line = line.split()
	Lx,Lz = float(line[0]),float(line[1])

print(Nx,Nz,Lx,Lz)

n_cores = 4

tempo_vec = []
topo_vec = []

for cont in range(0,10000,10):
    try:
        if n_cores>1:
            tempo = np.loadtxt("Tempo_"+str(cont)+".txt",unpack=True,delimiter=":",usecols=(1))
            tempo = tempo[0]
        else:
            tempo = np.loadtxt("Tempo_"+str(cont)+".txt",unpack=True,delimiter=":")
    except:
        print("Nao encontrou passo %d"%(cont))
        break
        

    x=[]
    z=[]
    c=[]
    s=[]
    for rank in range(n_cores):
        x1,z1,c0,c1,c2 = np.loadtxt("step_"+str(cont)+"-rank_new"+str(rank)+".txt",unpack=True)

        cor =  (0,0,0)
        cor2 = (0,0,0)
        cor3 = (0,0,0)
        #print(cor)
		
        c = np.append(c,c0)
        x = np.append(x,x1)
        z = np.append(z,z1)
        s = np.append(s,c1)

    x = x[c>9999]
    z = z[c>9999]
    c = c[c>9999]
    s = s[c>9999]

    if cont==0:
        condition = (z>-155.0E3)&(z<-150.0E3)
    else:
        condition = (s==3)&(z>-170.0E3)&(z<-140.0E3)
    xx = x[condition]
    zz = z[condition]
    cc = c[condition]
    
    ic = np.argsort(cc)

    xx = xx[ic]
    zz = zz[ic]
    cc = cc[ic]
    
    if cont==0:
        xx0 = np.copy(xx)
        zz0 = np.copy(zz)
        cc0 = np.copy(cc)
        num_pontos = np.size(cc0)

    cc0_indices = np.where(np.isin(cc0,cc)) 
    cc_indices = np.where(np.isin(cc,cc0)) 

    dx = xx[cc_indices]-xx0[cc0_indices]
    dz = zz[cc_indices]-zz0[cc0_indices]

    tempo_vec = np.append(tempo_vec,tempo)
    topo_vec = np.append(topo_vec,np.max(dz))

    print(cont,np.size(xx))

plt.close()

#plt.plot(xx[cc_indices],dz,".")
plt.plot(tempo_vec/1.0E6,topo_vec,".r",label="MANDYOC")


#plt.title("%.2f"%(tempo/1.0E6))

#plt.ylim([-300,1000])
plt.grid("on")

x,y = np.loadtxt("UNDERWORLD.txt",unpack=True, delimiter = ",")
plt.plot(x,y,"g",label="UNDERWORLD")

x,y = np.loadtxt("STAGYY.txt",unpack=True, delimiter = ",")
plt.plot(x,y,"k",label="STAGYY")

x,y = np.loadtxt("I2VES.txt",unpack=True, delimiter = ",")
plt.plot(x,y,"b",label="I2VIS")

plt.legend(loc=4)

plt.ylim(0,900)
plt.xlim(0,20)

plt.savefig("cmax_topo.png")    

plt.xlim(0,0.8)
plt.ylim(0,400)

plt.legend(loc=4)

plt.savefig("cmax_topo_zoom.png")    
    
