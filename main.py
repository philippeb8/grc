from __future__ import division, print_function
import numpy as np
from numpy.random import randn
from numpy.fft import rfft
from numpy import asarray
from scipy import signal
from scipy.misc import derivative
from PIL import Image
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

galaxies = {}

file1 = open('db/0/gm.txt', 'r')
  
while True:
    line = file1.readline()
  
    if not line:
        break
    
    fields = line.split()
    
    galaxies[fields[0]] = float(fields[9]) * 1e9 * float(fields[1]) * 2e30
  
file1.close()

velocities = {}

file2 = open('db/0/grc.txt', 'r')
  
while True:
    line = file2.readline()
  
    if not line:
        break
    
    fields = line.split()
    
    if fields[0] not in velocities:
        velocities[fields[0]] = []
    
    velocities[fields[0]].append((float(fields[2])*3.08567758128e+19, float(fields[3])*1000))
  
file2.close()

plt.subplot(1, 1, 1)
plt.grid(True, which='both')

mng = plt.get_current_fig_manager()
mng.resize(1024, 768*2/3)

G = 6.674e-11

totalstddev = 0.0
sortedgalaxies = {}

for count, key in enumerate(galaxies):
    print(math.floor(count * 100 / len(galaxies)), "%", end='\r')
    sig_ff = []

    graph = velocities[key]
    
    xvelocities = [x for (x, y) in graph]
    yvelocities = [y for (x, y) in graph]
    
    shortestm = 0.0
    shortestw = 0.0
    shortesth = 0.0
    shorteststddev = 9e99
    shortestvelocities = []
    
    for h in np.arange(1e20 / 100, 1e20 * 100, 1e20 / 2):
        for m in np.arange(galaxies[key] / 15, galaxies[key] * 15, pow(10, math.floor(math.log(galaxies[key], 10))) / 2):
            ftvelocities = [h/(m/x+h)*math.sqrt(G*m/x) for x in xvelocities]
            w = (yvelocities[-1] - ftvelocities[-1]) / xvelocities[-1]
            ft2velocities = [y + w * x for (x, y) in zip(xvelocities, ftvelocities)]
            
            stddev = np.std([(y1 - y2) for (y1, y2) in zip(yvelocities, ft2velocities)])
            
            if stddev < shorteststddev:
                shortestm = m
                shortestw = w
                shortesth = h
                shorteststddev = stddev
                shortestvelocities = ft2velocities
   
    totalstddev += shorteststddev
    sortedgalaxies[shorteststddev] = (key, galaxies[key], shortestm, shortestw, shortesth, xvelocities, yvelocities, shortestvelocities)

print("Mean stddev =", totalstddev / len(galaxies), "m/s")

for count, (key, value) in enumerate(sorted(sortedgalaxies.items(), key=lambda x: x[0])):
    observed = [value[5], value[6]]
    theoretical = [value[5], value[7]]

    plt.clf()
    plt.title(str(count + 1) + "/" + str(len(sortedgalaxies)) + ": " + str(value[0]) + " Tangential Velocity vs Radius (h = " + '{:.2e}'.format(value[4]) + " kg/m, m/m_0 = " + '{:.2e}'.format(value[2] / value[1]) + ")")

    plt.plot(observed[0], observed[1], label='Observed Data')
    plt.plot(theoretical[0], theoretical[1], label='Theoretical Data')
    plt.legend(loc="upper center")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    if not plt.waitforbuttonpress():
        exit()
