#/usr/bin/python -env

import numpy as np
import matplotlib.pyplot as plt
from rdf import pairCorrelationFunction_3D
import csv

# Read input file
with open("inputRdf", "rb") as f:
    reader = csv.reader(f, delimiter="\t")
    pars = []
    for row in reader:
        par = row[1]
        pars.append(par)     	
	
# Read particle data
dat = np.loadtxt(pars[0])
x = dat[:,0]
y = dat[:,1]
z = dat[:,2]
num_particles = len(x)
print "Number of particles =",num_particles

# Particles diameter
particle_radius = float(pars[1])/2.

# Bin size
dr = particle_radius/float(pars[2])

# Maximum radius
rMax = float(pars[3])*particle_radius

# Domain size
Lx = float(pars[4])
Ly = float(pars[5])
Lz = float(pars[6])

# Compute pair correlation
g_r, r, reference_indices = pairCorrelationFunction_3D(x, y, z, Lx, Ly, Lz, rMax, dr)

# Write into file
f = open('RDF.dat', 'w')
for i in xrange(len(g_r)):
    f.write(str(r[i]/(2.*particle_radius)) + " " + str(g_r[i]) + "\n")
f.close()

# Plot
plt.figure()
plt.plot(r/(2.*particle_radius), g_r, color='black')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.xlim( (0, rMax/(2.*particle_radius)) )
plt.ylim( (0, 1.05 * g_r.max()) )
plt.show()


