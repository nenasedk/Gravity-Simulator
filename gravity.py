from oct_tree import *
import numpy as np
from numpy.linalg import norm
from scipy import special
import os 
import sys

class Gravity:
    def __init__(self):
        return

    def read_in_data(self,filepath, filename):
        if not filepath.endswith('/'):
            filepath += '/'
        if not filename.endswith('.txt') or not filename.endswith('.ascii'):
            print("Can only read in .txt or .ascii data files.")
        
        with open(filepath + filename) as f:
            header = list(map(int,f.readline().split()))
            f.close()
        data = np.genfromtxt(filepath + filename,skip_header=1)
        data = data.reshape(9,header[0])
        return header, data

    def density(self,mass,posn,n):
        coords = self.cart_to_sphere(posn)
        coordinds = coords[:,0].argsort()
        s_coords = coords[:,0][coordinds]
        
        s_mass = mass[coordinds]
        max_rad = np.max(s_coords)
        h = max_rad/n
        radius = 0.0
        radout = np.zeros(n-1)
        density = np.zeros(n-1)
        j = 0
        for i in range(n-1):
            volume = (4.0/3.0)*np.pi * (((i+1)*h) - (i*h))**3
            m_shell = 0.0
            while(s_coords[j]<(i+1)*h):
                m_shell += s_mass[j]
                j+=1
            density[i] = m_shell/volume
            radout[i] = i*h + (i+1)*h/2.0
        return radout,density
            
    def analytical_force(self):
        radius = np.arange(0.0001,1.0,0.0001)
        return

    def relaxation_time(self,nParticles,t_cross):
        t_relax = nParticles * t_cross/(8.0*np.log(nParticles))
        return t_relax

    def crossing_time(self,mass,posn,vels):
        M_t = np.sum(mass) # Compute total mass
        coords = self.cart_to_sphere(posn)
        # Sort by radius
        coordinds = coords[:,0].argsort() 
        s_coords = coords[:,0][coordinds[::0]]
        s_mass = mass[coordinds[::0]]
        m = 0.0
        i = 0 
        while(m<M_t/2): # Find Half-Mass Radius
            m += s_mass[i]
            i += 1            
        r_half = s_coords[i]
        V_p = np.sqrt(M_t/(2.0 * r_half)) # Velocity
        return r_half/v_p # Return crossing time.
    
    def cart_to_sphere(self,data):
        newcoords = np.zeros(data.shape, dtype  = np.float64)
        dist = np.square(data[:,0]) + np.square(data[:,1]) + np.square(data[:,2])
        newcoords[:,0] = np.sqrt(np.square(data[:,0]) + np.square(data[:,1]) + np.square(data[:,2]))
        newcoords[:,1] = np.arctan2(data[:,1],data[:,0])
        newcoords[:,2] = np.arccos(data[:,2]/newcoords[:,0])
        return newcoords  
    

    
class DirectNBody(Gravity):
    def __init__(self,
                 filepath = None,
                 filename = "data.ascii"):
        if filepath is None:
            self.filepath = os.getcwd()
        else:
            self.filepath = filepath
        self.filename = filename
        # NOTE Read in from file indexed from 1, not 0
        self.nParticles = -1
        self.nGasParticles = -1
        self.nStarParticles = -1

        self.mass = np.zeros(0)
        self.posn = np.zeros((3,0))
        self.vels = np.zeros((3,0))
        self.soft = np.zeros(0)
        self.pot = np.zeros(0)
        self.s_soft = 0.0

    """
    Utility Functions
    """
    def set_nParticles(self,n):
        self.nParticles = n
    def set_nGasParticles(self,ng):
        self.nGasParticles = ng
    def set_nStarParticles(self,ns):
        self.nStarParticles = ns
    def set_mass(self,mass):
        self.mass = mass
    def set_posn(self,ps):
        self.posn = ps
    def set_vels(self,vs):
        self.vels = vs
    def set_soft(self,s):
        self.soft = s
    def set_pot(self,p):
        self.pot = p
        
    def set_variables(self):
        header,data = self.read_in_data(self.filepath,self.filename)
        self.nParticles = header[0]
        self.nGasParticles = header[1]
        self.nStarParticles = header[2]

        self.mass = np.zeros(self.nParticles)
        self.posn = np.zeros((3,self.nParticles))
        self.vels = np.zeros((3,self.nParticles))
        self.soft = np.zeros(self.nParticles)
        self.pot = np.zeros(self.nParticles)
        
        self.mass = data[0]
        self.posn[0] = data[1]
        self.posn[1] = data[2]
        self.posn[2] = data[3]
        self.vels[0] = data[4]
        self.vels[1] = data[5]
        self.vels[2] = data[6]
        self.soft = data[7]
        self.pot = data[8]

        self.posn = self.posn.reshape((self.nParticles,3))
        self.vels = self.vels.reshape((self.nParticles,3))
        return
    
    def direct_force(self,mass,posn,soft):
        try:
            assert(self.nParticles > 0)
        except AssertionError:
            print("Please input data, no particles to compute!")
            sys.exit(1)
        try:
            assert(len(mass) == self.nParticles)
        except AssertionError:
            print("Incorrect number of masses")
            sys.exit(1)
        try:
            assert(len(posn) == self.nParticles)
        except AssertionError:
            print("Incorrect number of positions, " + str(len(posn)) + " positions for " + str(self.nParticles) + " particles.")
            sys.exit(1)

        force = np.zeros((self.nParticles ,3))
        #m2 = np.outer(mass,mass)
        if(np.isscalar(soft)):
            # Naive implementation, should be able to stop at nparticles/2?
            for i in range(self.nParticles):
                for j in range(self.nParticles):
                    if i == j:
                        continue
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j]) + soft)**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
        else:
            for i in range(self.nParticles-1):
                for j in range(self.nParticles-1):
                    if i == j:
                        continue
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j]) + soft[i])**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
                    
        return force

    def direct_force_fast(self,mass,posn,soft):
        try:
            assert(self.nParticles > 0)
        except AssertionError:
            print("Please input data, no particles to compute!")
            sys.exit(1)
        try:
            assert(len(mass) == self.nParticles)
        except AssertionError:
            print("Incorrect number of masses")
            sys.exit(1)
        try:
            assert(len(posn) == self.nParticles)
        except AssertionError:
            print("Incorrect number of positions, " + str(len(posn)) + " positions for " + str(self.nParticles) + " particles.")
            sys.exit(1)

        force = np.zeros((self.nParticles ,3))
        print(len(mass),len(posn))
        # m2 = np.outer(mass,mass)
        if(np.isscalar(soft)):
            # Naive implementation, should be able to stop at nparticles/2?
            for i in range(0,5):
                if(i%500 == 0):
                    print(str(i/500) + "%")
                for j in range(self.nParticles-1,i,-1):
                    if i == j:
                        break
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j]) + soft)**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    force[j] = force[j] - force[i]
        else:
            for i in range(0,self.nParticles-1):
                if(i%500 == 0):
                    print(str(i/500) + "%")
                for j in range(self.nParticles-1,i,-1):
                    if i == j:
                        break
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j]) + soft[i])**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    force[j] = force[j] - force[i]
                    
        return force


    
class OctGrid(Gravity):
    def __init__(self,
                 filepath = None,
                 filename = "data.ascii"):
        if filepath is None:
            self.filepath = os.getcwd()
        else:
            self.filepath = filepath
        self.filename = filename
        # NOTE Read in from file indexed from 1, not 0
        self.nParticles = -1
        self.nGasParticles = -1
        self.nStarParticles = -1

        self.tree = None
        self.size = 0
        self.soft = np.zeros(0)
        self.pot = np.zeros(0)
        self.s_soft = 0.0

    """
    Utility Functions
    """
    def set_nParticles(self,n):
        self.nParticles = n
    def set_nGasParticles(self,ng):
        self.nGasParticles = ng
    def set_nStarParticles(self,ns):
        self.nStarParticles = ns
    def set_soft(self,s):
        self.s_soft = s
    def set_pot(self,p):
        self.pot = p
    def set_size(self,s):
        self.size = s
        
    def buildParticleList(self):
        header,data = self.read_in_data(self.filepath,self.filename)
        self.nParticles = int(header[0])
        self.nGasParticles = int(header[1])
        self.nStarParticles = int(header[2])
        particles = []
        self.size = np.max([np.max(data[1]),
                            np.max(data[2]),
                            np.max(data[3])])*2.0
        for i in range(self.nParticles):
            particles.append(Particle(np.array([data[1][i],data[2][i],data[3][i]]),
                                      np.array([data[4][i],data[5][i],data[6][i]]),
                                      data[0][i]))
        self.soft = data[7]
        self.pot = data[8]
        return particles

    def buildTree(self,particles):
        self.tree = OctTree(self.size,1)
        self.tree.buildTree(particles)
        return
    
    def close_enough(p1,p2):
        return

    def multipole(self, p1, p2tree, L):
        r1 = p1.posn
        CoM = p2tree.GetCenterOfMass()
        MT = p2tree.GetTotalMass()
        s = r - CoM #Maybe need the actual vect?
        s_sph = CartToSphere(s)
        # cos_theta = np.dot(r1,r2)/(norm(r1) * norm(r2))
        # phi = (1/norm(r2))*np.sum([((norm(r1)/norm(r2))**n)*special.legendre(n)(cos_theta) for n in range(N)])
        multipoles = []
        moments = self.calcMultipoleMoments(p1,p2tree,L,M)
        for l in range(L+1):
            phi_l = 0
            for m in range(-l,l+1):
                sph = sph_harm(m,l,s_sph[1],s_sph[2])
                Q = moments[(l,m)]
                phi_l = np.sqrt(4*np.pi/(2*l+1)) * Q * sph /(s_sph[0]**(l+1))
            multipoles.append(phi_l)
        pot = sum(multipoles[:L+1])
        return pot

    def shiftMultipoleMoments(self,moments):
        return
    
    def calcMultipoleMoments(self,p1,p2tree,L,M):
        moments = {}
        for l in range(L+1):
            phi_l = 0
            for m in range(-l,l+1):
                moments[(l,m)] = self.calcMultipoleCoefficient(p2tree,l,m)        
        return moments
    
    def calcMultipoleCoefficient(self,node,l,m):
        Q = 0
        for p in node.data:
            r2 = CartToSphere(p.posn)
            sph = special.sph_harm(m,l,r2[1],r2[2])
            Q += p.mass * r2[0] * np.conj(sph)
        return Q.real
    def potential():
        return
    
    def compute_force():
        return
    
