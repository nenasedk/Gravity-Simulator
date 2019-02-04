from gravity import *
import numpy as np
import matplotlib.pyplot as plt
import time as t
#import os
#import sys
"""
soft_scale = 100.0
def plummerdens(d0,r0,r):
    d = d0/(1 + (r**2 / (20 * r0**2)))
    return d

def Hernquist(M,a,r):
    d = M/(2*np.pi) *(a/r) * 1/(r + a)**3
    return d

nBody = DirectNBody()
nBody.set_variables()
nBody.set_nParticles(1000)
polar = nBody.cart_to_sphere(nBody.posn)
eps = np.std(polar[:,0])

t_start = t.time()
#force = nBody.direct_force(nBody.mass,nBody.posn,eps)
force2 = nBody.direct_force_fast(nBody.mass,nBody.posn,eps*soft_scale)
t_end = t.time()


myFile = open("small_force_" + str(nBody.nParticles) + "_" + str(soft_scale) + ".txt",'w+')
myFile.write(str(nBody.nParticles) + "," + str(t_end - t_start) + '\n')
for i in range(0,len(force2)):
    ostring = str(force2[i][0]) + "," + str(force2[i][1]) + "," + str(force2[i][2]) + "," +str(polar[i][0]) + "," + str(nBody.posn[i][0]) + ","+ str(nBody.posn[i][1]) + ","+ str(nBody.posn[i][2]) + '\n'
    myFile.write(ostring)
myFile.close()
"""

fig,ax = plt.subplots(figsize = (12,8))
filedir = "1sigma/"
time = []
n = []
for filename in os.listdir(filedir):
    with open(filedir + filename) as f:
        header = list(map(float,f.readline().split(',')))
        f.close()
    #sp = filename.split('_')
    #sp2 = sp[3].split('t')
    data = np.genfromtxt(filedir+filename,delimiter = ',',skip_header =1)
   
    force = np.array([data[:,0],data[:,1],data[:,2]]).T
    n.append(header[0])
    time.append(header[1])
        
    radius = data[:,3]
    #f = -1*np.linalg.norm(force2,axis = 1)
    #ax.scatter(radius,np.linalg.norm(force,axis = 1),marker = '.', s = 2,label = str(len(force)) + " Particles",linewidth = 3)
    #ax.scatter(radius,np.linalg.norm(force,axis = 1),marker = '.', s = 2,label =  sp2[0][:-1],linewidth = 3)
y = [x for _,x in sorted(zip(n,time))]

fit = np.polyfit(np.array(sorted(n)),np.array(y),2)
fitfn = np.poly1d(fit)
ti = np.linspace(10.0,max(n),500000)
ax.scatter(sorted(n),y,marker = '.',s = 25,linewidth = 2, c = u'#ff7f0e')
ax.plot(ti,fitfn(ti))
#ax.set_xlim(0,40)
ax.set_yscale("log")
ax.set_ylim(1e-4,5e5)
ax.set_title("O(N$^{2}$) runtime for direct n-body force calculation with quadratic fit")
ax.set_xlabel("Number of Particles")
ax.set_ylabel("Run time [s]")
#ax.legend(title = "Softening",loc = 'upper right')
plt.show()

'''
lin_rad,lin_density,lin_err,s_coords = nBody.density(nBody.mass,nBody.posn,5001,False)
log_rad,log_density,log_err,a = nBody.density(nBody.mass,nBody.posn,25,True)
h_rad,h_density,h_err,b = nBody.density(nBody.mass,nBody.posn,5000,False,True)

print(np.min(log_rad),np.min(lin_rad),np.min(h_rad),np.min(s_coords))
mr = np.max(s_coords)
minr = np.min(s_coords)
r = np.linspace(minr/mr,0.01,1000000)

s_lin = minr/np.min(lin_rad)
s_log = minr/np.min(log_rad)
s_h = minr/np.min(h_rad)
dp = plummerdens(1.0,minr/mr,r)
d = Hernquist(1.0,mr/(np.pi),r)

fig,ax = plt.subplots(figsize = (12,8))
ax.plot(r,d,'--',label = "Hernquist Model",linewidth = 3)
ax.plot(r,dp,'--', label = "Plummer-like Model",linewidth = 3)

ax.errorbar(s_log * log_rad/mr,
            log_density/np.max(log_density),
            yerr = log_density/(np.max(log_density)*log_err),
            marker = '.',
            linewidth = 1,
            markersize = 4,
            alpha = 0.8,
            label = "Log-binned data")
ax.errorbar(s_lin*lin_rad/mr,
            lin_density/np.max(lin_density),
            yerr = lin_density/(np.max(lin_density)*lin_err),
            marker = 'o',
            markersize = 2,
            linewidth = 1,
            alpha = 0.8,
            label = "Linear-binned data")
ax.errorbar(s_h*h_rad/mr,
            h_density/np.max(h_density),
            yerr = h_density/(np.max(h_density)*h_err),
            marker = 's',
            markersize = 4,
            linewidth = 1,
            alpha = 0.8,
            label = "Even-binned data")

legend = ax.legend(loc = 'upper right')
ax.set_title("Radial density profiles of a spherical galaxy.",fontsize = 24)
ax.set_xlabel("Normalised Radius",fontsize = 18)
ax.set_ylabel("Normalised Density",fontsize = 18)
ax.set_xscale("log")
plt.show()
'''
