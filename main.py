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

file1 = open('gm.txt', 'r')
  
while True:
    line = file1.readline()
  
    if not line:
        break
    
    fields = line.split()
    
    galaxies[fields[0]] = float(fields[9]) * 1e9 * float(fields[1]) * 2e30
  
file1.close()

velocities = {}

file2 = open('grc.txt', 'r')
  
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

plt.show(False)

G = 6.674e-11

for key in galaxies:
    sig_ff = []

    graph = velocities[key]
    
    xvelocities = [x for (x, y) in graph]
    yvelocities = [y for (x, y) in graph]
    
    # Skipping low speed galaxies
    if max(yvelocities) < 100000:
        continue

    sig_ff.append([xvelocities, yvelocities])
    
    shortestm = 0.0
    shortestw = 0.0
    shortesth = 0.0
    shortestdiff = 9e99
    shortestvelocities = []

    for m in np.arange(galaxies[key] / 10, galaxies[key] * 10, pow(10, math.floor(math.log(galaxies[key], 10)))):
        for h in np.arange(1e19, 5e20, 1e19):
            ftvelocities = [h/(m/x+h)*math.sqrt(G*m/x) for x in xvelocities]
            w = (yvelocities[-1] - ftvelocities[-1]) / xvelocities[-1]
            ft2velocities = [y + w * x for (x, y) in zip(xvelocities, ftvelocities)]
            
            diff = sum([abs(y1 - y2) for (y1, y2) in zip(yvelocities, ft2velocities)])
            
            if diff < shortestdiff:
                shortestm = m
                shortestw = w
                shortesth = h
                shortestdiff = diff
                shortestvelocities = ft2velocities

    sig_ff.append([xvelocities, shortestvelocities])

    plt.clf()
    plt.title(key + " Rotation Curve vs Radius (h = " + str(shortesth) + ", m/m_0 = " + str(shortestm / galaxies[key]) + ")")

    plt.plot(sig_ff[0][0], sig_ff[0][1], label='Observed Data')
    plt.plot(sig_ff[1][0], sig_ff[1][1], label='Theoretical Data')
    plt.legend(loc="upper center")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    if not plt.waitforbuttonpress():
        exit()
