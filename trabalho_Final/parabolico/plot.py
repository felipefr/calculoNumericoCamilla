import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("saida.txt",skiprows=2)

f = open('saida.txt')

Ntempos = int(f.readline())
N = int(f.readline())

data = data.reshape((Ntempos,N))

Nx = 21
Ny = 21

x = np.linspace(0.0,1.0,Nx)
y = np.linspace(0.0,1.0,Ny)

X, Y = np.meshgrid(x,y)

plt.figure(1)

for n in range(Ntempos): 
#~ plt.subplot(121)
	
	Un = data[n,:].reshape((Nx,Ny))
	plt.imshow(U,interpolation = "bilinear", cmap = 'jet')
	#~ plt.xticks(x)
	plt.colorbar()
	plt.grid()

	plt.savefig('placa_' + str(n) + '.png')

#~ plt.subplot(122)
#~ plt.imshow(Uex,interpolation = "bilinear", cmap = 'jet')
#~ plt.colorbar()
#~ plt.grid()

#~ plt.show()

