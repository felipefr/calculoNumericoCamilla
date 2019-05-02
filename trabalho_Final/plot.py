import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("saida.txt")

Nx = 20
Ny = 20
x = np.linspace(0.0,1.0,Nx)
y = np.linspace(0.0,1.0,Ny)

X, Y = np.meshgrid(x,y)

U = data.reshape((Nx,Ny))
Uex = np.cos(X) + np.sin(Y)

plt.figure(1)
plt.subplot(121)
plt.imshow(U,interpolation = "bilinear", cmap = 'jet')
plt.xticks(x)
plt.colorbar()
plt.grid()

plt.subplot(122)
plt.imshow(Uex,interpolation = "bilinear", cmap = 'jet')
plt.colorbar()
plt.grid()

plt.show()

