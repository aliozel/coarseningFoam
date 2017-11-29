import numpy as np
import matplotlib.pyplot as plt
#from utilities import *
from rdf import pairCorrelationFunction_3D

# Particle setup
Lx = 20.0
Ly = 20.0
Lz = 20.0
num_particles = 10000

# Calculation setup
dr = 0.1

### Random arrangement of particles ###
particle_radius = 0.1
rMax = 5*particle_radius
x = np.random.uniform(low=0, high=Lx, size=num_particles)
y = np.random.uniform(low=0, high=Ly, size=num_particles)
z = np.random.uniform(low=0, high=Lz, size=num_particles)

# Compute pair correlation
g_r, r, reference_indices = pairCorrelationFunction_3D(x, y, z, Lx, Ly, Lz, rMax, dr)

# Visualize
plt.figure()
plt.plot(r, g_r, color='black')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.xlim( (0, rMax) )
plt.ylim( (0, 1.05 * g_r.max()) )
plt.show()