import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import glob


if len(sys.argv) > 1:
	step_initial = int(sys.argv[1])
else:
	step_initial = 0

if len(sys.argv) > 2:
	step_final = int(sys.argv[2])
else:
	step_final = 2100

if len(sys.argv) > 3:
	d_step = int(sys.argv[3])
else:
	d_step = 10

with open("param.txt","r") as f:
	line = f.readline()
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


for cont in range(step_initial,step_final,d_step):#
	print(cont)

	A = np.loadtxt("time_"+str(cont)+".txt",dtype='str')
	AA = A[:,2:]
	AAA = AA.astype("float")
	time = np.copy(AAA)[0, 0] / 1.0e6

	print(f'Time = {time:.1f} Myr\n\n')

	A = pd.read_csv("density_"+str(cont)+".txt",delimiter = " ",comment="P",skiprows=2,header=None)
	A = A.to_numpy()
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Nz),order='F')
	TTT = TT[:,:]
	TTT = np.transpose(TTT)
	rho = np.copy(TTT)


	A = pd.read_csv("viscosity_"+str(cont)+".txt",delimiter = " ",comment="P",skiprows=2,header=None)
	A = A.to_numpy()
	TT = A*1.0
	TT = np.reshape(TT,(Nx,Nz),order='F')
	TTT = np.transpose(TT)
	TTT = np.log10(TTT)

	viscosity = np.copy(TTT)
	viscosity[viscosity==-np.inf] = 0


	plt.figure(figsize=(20,5))

	cr = 255.
	color_plume = (207./cr,226./cr,205./cr)
	color_mantle = (155./cr,194./cr,155./cr)
	color_lid = (228./cr,156./cr,124./cr)

	plt.contourf(xx, zz, viscosity,
		levels=[19, 20, 21, 23],
		colors=[color_plume, color_mantle, color_lid]
	)

	plt.colorbar()
	plt.title(f'Time = {time:.1f} Myr\n\n')

	plt.savefig("fig_viscosity_{:05}.png".format(cont*1))
	plt.close()
