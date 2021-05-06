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

xi = np.linspace(0,Lx/1000,Nx)
zi = np.linspace(-Lz/1000,0,Nz)
xx,zz = np.meshgrid(xi,zi)

n_cores = 4

for cont in range(0,4000,40):
	print(cont)
	try:
		if n_cores>1:
			tempo = np.loadtxt("time_"+str(cont)+".txt",unpack=True,delimiter=":",usecols=(1))
			tempo = tempo[0]
		else:
			tempo = np.loadtxt("time_"+str(cont)+".txt",unpack=True,delimiter=":")
	except:
		print("Step %d not found"%(cont))
		break
	
	A = np.loadtxt("density_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Nz),order='F')
	TTT = TT[:,:]
	plt.close()
	plt.figure(figsize=(2,2))
	plt.imshow(np.transpose(TTT[:,::-1]), extent = [0,0.9142,0,1.0])
	plt.title("t' = %.1f"%(tempo*365*24*3600.*1.0E-12)) #adimensional time
	#tempo in years
	# t' = tempo * seconds_per_year * kappa / h^2
	# t' = tempo * seconds_per_year * (rho g h^3/ Rb visc) / h^2
	
	plt.savefig("density_{:05}.png".format(cont))
	
	




