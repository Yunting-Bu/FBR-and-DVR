import numpy as np
import matplotlib.pyplot as plt
import math

a = -10.0
b = 10.0

def sin_basis(x,n):
    phi = math.sqrt(2/(b-a))*math.sin(n*math.pi*(x-a)/(b-a))
    return phi

ngrid = int(input('The number of grids:'))

xl = np.array([a + l*(b-a)/(ngrid+1) for l in range(1, ngrid+1)])

i = int(input(f'level(the value of n, from 0 to {ngrid-1}): '))

FBR = 'FBR_coeff.log'
DVR = 'DVR_coeff.log'
FBRcoeff = np.loadtxt(FBR)
DVRcoeff = np.loadtxt(DVR)

x_FBR = np.linspace(-10, 10, 200)
y_FBR = np.array([
    sum(FBRcoeff[n, i] * sin_basis(x, n + 1) for n in range(ngrid))
    for x in x_FBR
])

x_DVR = xl
y_DVR = np.array([DVRcoeff[l, i] for l in range(ngrid)])

plt.figure(figsize=(10, 6))
plt.plot(x_FBR, y_FBR, color='blue')
plt.xlabel('x')
plt.ylabel('Psi(x)')
plt.title(f'FBR wavefunction for level {i+1}')
plt.savefig('FBRwavefunction.png')

plt.figure(figsize=(10, 6))
plt.plot(x_DVR, y_DVR, color='red') 
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'DVR wavefunction for level {i+1}')
plt.savefig('DVRwavefunction.png')