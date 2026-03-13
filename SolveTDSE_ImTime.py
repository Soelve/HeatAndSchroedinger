"""
 This script estimates the ground state of a one dimensional quantum
 system by solving the time-dependent Schrödinger equation in imaginary
 time. The Hamiltonian is estimated using a three point approximation for 
 the kinetic energy operator.

 Numerical inputs:
   L       - The extension of the spatial grid 
   N       - The number of grid points
   Tfinal  - The duration of the simulation
   dt      - The step size in time

 Physical inputs:
   V0      - The height of the barrier (should be negative)
   w       - The width of the barrier.
   s       - A smoothness parameter for the barrier.
 
 There is also an input for adjusting the pause between each frame displayed
 during the simulation.  
 All inputs are hard coded initially.
"""   

# Import libraries  
import numpy as np
from matplotlib import pyplot as plt

# Numerical grid parameters
L = 100
N = 250

# Numerical time parameters
Tfinal = 15
dt = 0.1

# Pause - used to adjust the speed of the "animation"
PauseVal = 0.02

# Inputs for the potential
V0 = -1
w = 4
s = 5

# Define function for the barrier
def Vpot(x):
    return V0/(np.exp(s*(np.abs(x)-w/2))+1)


# Set up grid
x = np.linspace(-L/2, L/2, N+1)
dx = L/N

# Set up matrix for the double derivative - 3 point formula
# Allocate and declare
H = np.zeros((N+1, N+1), dtype=complex)
# Endpoints
H[0, 0:2] = [-2, 1]
H[N, (N-1):(N+1)] = [1, -2]
# Interior points
for n in range (1,N):
  H[n, [n-1, n, n+1]] = [1, -2, 1]
# Correct pre-factors
H = H/dx**2
H = -1/2*H

# Add potential
H = H + np.diag(Vpot(x))

# Crank-Nicolson propagator
Imat = np.eye(N+1)
U_CN = np.matmul(np.linalg.inv(Imat + dt/2*H), Imat - dt/2*H)

# Set up initial Gaussian - analytically
Psi0 = np.random.rand(N+1)
Norm = dx*np.sum(np.abs(Psi0)**2)
Psi0 = Psi0/np.sqrt(Norm)

# Initiate plots
plt.ion()
fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot()
# Initiate plot of wave function
line1, = ax.plot(x, np.abs(Psi0)**2, '-', color='black')
line2, = ax.plot(x, 0.1*Vpot(x)/np.abs(V0), '-', color='red')
ax.set(ylim=(-.1, .5))
plt.grid(visible=True)


# Initiate wave functons and time
# Turn Psi0 into (N+1) \times 1 matrix (column vector)
Psi = np.matrix(Psi0)                 
Psi = Psi.reshape(N+1,1)
t = 0

# Vector with energy estimates
tVector = np.arange(0, Tfinal, dt)
EnergyVector = np.zeros(len(tVector))
index = 0

# Iterate while time t is less than final time Tfinal
for t in tVector:
  # Update numerical wave function
  Psi = np.matmul(U_CN, Psi)
  
  # Renormalize
  Norm = dx*np.sum(np.power(np.abs(Psi),2))
  Psi = Psi/np.sqrt(Norm)

  # Update data for plots
  line1.set_ydata(np.power(np.abs(Psi),2))
  
  # Update plots
  fig.canvas.draw()
  #fig.canvas.flush_events()
  plt.pause(PauseVal)
    
  # Estimate energy
  E0 = -np.log(Norm)/(2*dt)
  EnergyVector[index] = E0
  index = index+1

# Plot energy estimates in "time"
plt.figure(2)
plt.clf()
plt.plot(tVector, EnergyVector, 'k-')
plt.grid(visible=True)
plt.show()

# Write final energy estimate
print(f'Ground state energy: {E0:.2f}')  