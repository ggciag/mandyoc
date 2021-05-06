import numpy as np
import matplotlib.pyplot as plt

Nx = 81
Nz = 81

Lx = 91.42
Lz = 100.0

xn = np.linspace(0,Lx,Nx)

Liton = -0.8*Lz + 0.02*Lz*np.cos(np.pi*xn/Lx)

Liton = Liton*1000

f = open("interfaces_creep.txt","w")

f.write("C 1.0 1.0\n")
f.write("rho -1000. 0.\n")
f.write("H 0.0E-12 0.0E-12\n")
f.write("A 0.0 0.0\n")
f.write("n 0.0 0.0\n")
f.write("Q 0.0 0.0\n")
f.write("V 0.0 0.0\n")

for i in np.arange(Nx):
	f.write("%lf\n"%(Liton[i]))

f.close()




