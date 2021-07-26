import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import glob


step_initial = int(sys.argv[1])
step_final = int(sys.argv[2])

if (len(sys.argv)>3): d_step = int(sys.argv[3])
else: d_step = 10

with open("param.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx = int(line[2])
	line = f.readline()
	line = line.split()
	Nz = int(line[2])
	line = f.readline()
	line = line.split()
	Lx = float(line[2])
	line = f.readline()
	line = line.split()
	Lz = float(line[2])

print(Nx,Nz,Lx,Lz) 

xi = np.linspace(0,Lx/1000,Nx)
zi = np.linspace(-Lz/1000,0,Nz)
xx,zz = np.meshgrid(xi,zi)

h_air = 40.0


for cont in range(step_initial,step_final,d_step):#
	print(cont)


	A = np.loadtxt("time_"+str(cont)+".txt",dtype='str')  
	AA = A[:,2:]
	AAA = AA.astype("float") 
	tempo = np.copy(AAA)
	
	print("Time = %.1lf Myr\n\n"%(tempo[0]/1.0E6))


	A = pd.read_csv("density_"+str(cont)+".txt",delimiter = " ",comment="P",skiprows=2,header=None) 
	A = A.to_numpy()
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Nz),order='F')
	TTT = TT[:,:]
	TTT = np.transpose(TTT)
	rho = np.copy(TTT)


	A = pd.read_csv("strain_"+str(cont)+".txt",delimiter = " ",comment="P",skiprows=2,header=None) 
	A = A.to_numpy()
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Nz),order='F')
	TTT = np.transpose(TT)
	TTT[rho<200]=0
	TTT = np.log10(TTT)
	stc = np.copy(TTT)
	
	
	plt.close()
	plt.figure(figsize=(20,5))

	cr = 255.
	color_uc = (228./cr,156./cr,124./cr)
	color_lc = (240./cr,209./cr,188./cr)
	color_lit = (155./cr,194./cr,155./cr)
	color_ast = (207./cr,226./cr,205./cr)

	plt.contourf(xx,zz+h_air,rho,levels=[200.,2750,2900,3365,3900],
		colors=[color_uc,color_lc,color_lit,color_ast])


	print("stc",np.min(stc),np.max(stc))

	print("stc(log)",np.min(stc),np.max(stc))
	plt.imshow(stc[::-1,:],extent=[0,Lx/1000,-Lz/1000+h_air,h_air],
		zorder=100,alpha=0.2,cmap=plt.get_cmap("Greys"),vmin=-0.5,vmax=0.9)

	plt.text(100,10,"%.1lf Myr"%(tempo[0]/1.0E6))


	b1 = [0.74,0.41,0.2,0.2]
	bv1 = plt.axes(b1)

	A = np.zeros((100,10))

	A[:25,:]=2700
	A[25:50,:]=2800
	A[50:75,:]=3300
	A[75:100,:]=3400

	A = A[::-1,:]

	xA = np.linspace(-0.5,0.9,10)
	yA = np.linspace(0,1.5,100)

	xxA,yyA = np.meshgrid(xA,yA)
	air_threshold = 200
	plt.contourf(xxA,yyA,A,levels=[air_threshold,2750,2900,3365,3900],
			colors=[color_uc,color_lc,color_lit,color_ast])

	plt.imshow(xxA[::-1,:],extent=[-0.5,0.9,0,1.5],
			zorder=100,alpha=0.2,cmap=plt.get_cmap("Greys"),vmin=-0.5,vmax=0.9)

	bv1.set_yticklabels([])

	plt.xlabel(r"log$(\epsilon_{II})$",size=18)

	plt.savefig("litho_new_temper_{:05}.png".format(cont*1))
	




