from gravity import *
import numpy as np
import matplotlib.pyplot as plt
#import os
#import sys
def plummerdens(d0,r0,r):
    d = d0/(1 + (r**2 / (20 * r0**2)))
    return d

def Hernquist(M,a,r):
    d = M/(2*np.pi) *(a/r) * 1/(r + a)**3
    return d

nBody = DirectNBody()
nBody.set_variables()
nBody.set_nParticles(100)
polar = nBody.cart_to_sphere(nBody.posn)
eps = np.std(polar[:,0])

#force = nBody.direct_force(nBody.mass,nBody.posn,eps)
force2 = nBody.direct_force_fast(nBody.mass,nBody.posn,eps)

myFile = open("small_force.txt",'w+')
"""
for i in range(0,len(force2)):
    ostring = str(force2[i][0]) + "," + str(force2[i][1]) + "," + str(force2[i][2]) + "," +str(polar[i][0]) + "," + str(nBody.posn[i][0]) + ","+ str(nBody.posn[i][1]) + ","+ str(nBody.posn[i][2]) + '\n'
    myFile.write(ostring)
myFile.close()
"""
"""
data = np.genfromtxt("force.txt",delimiter = ',')

force = np.array([data[:,0],data[:,1],data[:,2]]).T
radius = data[:,3]
f = -1*np.linalg.norm(force,axis = 1)
"""
fig,ax = plt.subplots(figsize = (12,8))
ax.scatter(polar[:,0][0:len(force2)],np.linalg.norm(force2,axis = 1),marker = '.', s = 2,label = "Hernquist Model",linewidth = 3)
ax.set_xlim(0,40)
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
