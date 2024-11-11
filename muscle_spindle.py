import numpy as np
import os, math, operator, random, csv, scipy
import csv
from scipy.optimize import curve_fit
from scipy import special
from scipy import stats
from scipy.optimize import minimize
from pylab import *


h, simtime = 0.0001, 10 # both in seconds
n = np.int(simtime / h) + 1 # number of time steps


dx = np.hstack((np.ones(np.int(10 / h))*0.0, np.ones(np.int(5 / h)) * 1.0, np.ones(np.int(2 / h))*0.0))  # INPUT VELOCITY (mm/s). input here is time(s)/h for each velocity

T, U, R, X = np.ones(len(dx)), np.zeros(len(dx)), np.zeros(len(dx)), np.zeros(len(dx)) # time, spindle strain, spindle firing rate, muscle fiber length

a, b, c, g, p = 0.3, 250.0, -15.0, 350.0, 0.1 # primary deefferented constants
t0, u0, x0 = 0.0, 0.0374, -5.0 # time, strain in sensory zone, initial muscle length

#a, b, c, g, p = 50, 50.0, -20.0, 80.0, 0.1 # secondary deefferented constants
#t0, u0, x0 = 0.0, 0.305, -5.0 # time, strain in sensory zone, initial muscle length


for i in range(len(dx)):
    t1 = t0 + h # time
    x1 = x0 + h * (dx[i]) # position by integrating velocity
    u1 = u0 + h * (dx[i] - a * ((b * u0 - x0 + c)/(x0 - c - u0)) ** 3)
    r = g * (u0 + p * (dx[i] - a * ((b * u0 - x0 + c)/(x0 - c - u0)) ** 3))
    T[i], X[i], U[i], R[i]  = t1, x1, u1, r 
    t0, x0, u0, = t1, x1, u1
    print(i/len(dx)*100)


plot(T[95000::],R[95000::])
show()

plot(T[95000::],U[95000::])
show()

plot(T[95000::],X[95000::])
show()


# PLOT STATES OVER TIME
axiscolor = 'k'
titlefsize = 32
axisfsize = 28
scalefsize = 24 
fig, ax = plt.subplots(figsize=(16,8))
axis([-1,7,20,60])
plt.title('Muscle Spindle Firing Rate', fontsize = titlefsize)
plt.xlabel('Time (s)', fontsize = scalefsize)
plt.ylabel('r, (pps: peaks per second)',fontsize = scalefsize)
d1 = errorbar(T[95000::]-10.0, R[95000::], linestyle = '-', linewidth = 2.0, color = '#FD8B0B')
#lg = legend([d1], ['$\\theta (deg)$'], fontsize = 18, loc = 'lower center', ncol = 2)
#lg.get_frame().set_linewidth(0.0)
#lg.get_frame().set_alpha(0.0)
plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=False, direction='in')
plt.xticks(np.arange(-1.0, 8, 1.0))
plt.yticks(np.arange(20.0, 80, 20.0))
ax.spines['left'].set_color(axiscolor)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_color(axiscolor)
ax.tick_params(axis='x', colors=axiscolor, labelsize = scalefsize)
ax.tick_params(axis='y', colors=axiscolor, labelsize = scalefsize)
plt.tight_layout()
show()

# PLOT STATES OVER TIME
axiscolor = 'k'
titlefsize = 32
axisfsize = 28
scalefsize = 24 
fig, ax = plt.subplots(figsize=(16,8))
axis([-1,7,0.2,0.52])
plt.title('Nerve Ending Stretch', fontsize = titlefsize)
plt.xlabel('Time (s)', fontsize = scalefsize)
plt.ylabel('$\\mu (mm)$',fontsize = scalefsize)
d1 = errorbar(T[95000::]-10.0, U[95000::], linestyle = '-', linewidth = 2.0, color = '#FD8B0B')
#lg = legend([d1], ['$\\theta (deg)$'], fontsize = 18, loc = 'lower center', ncol = 2)
#lg.get_frame().set_linewidth(0.0)
#lg.get_frame().set_alpha(0.0)
plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=False, direction='in')
plt.xticks(np.arange(-1.0, 8, 1.0))
plt.yticks(np.arange(0.2, 0.6, 0.1))
ax.spines['left'].set_color(axiscolor)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_color(axiscolor)
ax.tick_params(axis='x', colors=axiscolor, labelsize = scalefsize)
ax.tick_params(axis='y', colors=axiscolor, labelsize = scalefsize)
plt.tight_layout()
show()


# PLOT STATES OVER TIME
axiscolor = 'k'
titlefsize = 32
axisfsize = 28
scalefsize = 24 
fig, ax = plt.subplots(figsize=(16,8))
axis([-1,7,-5.1,0.1])
plt.title('Muscle Length', fontsize = titlefsize)
plt.xlabel('Time (s)', fontsize = scalefsize)
plt.ylabel('x (mm)',fontsize = scalefsize)
d1 = errorbar(T[95000::]-10.0, X[95000::], linestyle = '-', linewidth = 2.0, color = '#FD8B0B')
#lg = legend([d1], ['$\\theta (deg)$'], fontsize = 18, loc = 'lower center', ncol = 2)
#lg.get_frame().set_linewidth(0.0)
#lg.get_frame().set_alpha(0.0)
plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False, bottom=True, top=False, left=True, right=False, direction='in')
plt.xticks(np.arange(-1.0, 8, 1.0))
plt.yticks(np.arange(-5.0, 1.0, 1.0))
ax.spines['left'].set_color(axiscolor)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_color(axiscolor)
ax.tick_params(axis='x', colors=axiscolor, labelsize = scalefsize)
ax.tick_params(axis='y', colors=axiscolor, labelsize = scalefsize)
plt.tight_layout()
show()


