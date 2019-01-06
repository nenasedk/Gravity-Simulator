from oct_tree import *
import numpy as np
from scipy import special
import os 
import sys

class Gravity:
    def __init__(self):
        return

    def read_in_data(self,filepath, filename):
        if not filepath.endswith('/'):
            filepath += '/'
        if not filename.endswith('.txt') or filename.endswith('.ascii'):
            print("Can only read in .txt or .ascii data files.")
        
        with open(filepath + filename) as f:
            header = list(map(int,f.readline().split()))
            f.close()
        data = np.genfromtxt(filepath + filename,skip_header=1)
        data.reshape(9,header[0] - 1)
        return header, data

    def density(self,mass,posn,n):
        coords = self.cart_to_sphere(posn)
        coordinds = coords[:,0].argsort()
        s_coords = coords[:,0][coordinds[::0]]
        s_mass = mass[coordinds[::0]]
        max_rad = np.max(s_coords)
        
        h = max_rad/n
        radius = 0.0
        density = np.zeros(n)
        j = 0
        for i in range(n-1):
            volume = (4.0/3.0)*np.pi * (((i+1)*h) - (i*h))**3
            m_shell = 0.0
            while(coords[j]<(i+1)*h):
                m_shell += s_mass[j]
                j+=1
            density[i] = m_shell/volume
        return density
            
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
        newcoords = np.array(data.shape)
        newcoords[:,0] = np.arctan2(data[:,1],data[:,0])
        newcoords[:,1] = np.sqrt(data[:,0]**2 + data[:,1]**2 + data[:,2]**2)
        newcoords[:,2] = np.arccos(data[:,2],newcoords[:,1])
        return newcoords  
    

    
class DirectNBody(Gravity):
    def __init__(self,
                 filepath = None,
                 filename = "data.ascii"):
        if filepath is None:
            self.filepath = os.cwd()
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
        
    def set_variables():
        header,data = self.read_in_data(self.filepath,self.filename)
        self.nParticles = header[0]
        self.nGasParticles = header[1]
        self.nStarParticles = header[2]
        
        self.mass = data[0]
        self.posn[0] = data[1]
        self.posn[1] = data[2]
        self.posn[2] = data[3]
        self.vels[0] = data[4]
        self.vels[1] = data[5]
        self.vels[2] = data[6]
        self.soft = data[7]
        self.pot = data[8]

        self.posn.reshape(self.nParticles - 1,3)
        self.vels.reshape(self.nParticles - 1 ,3)
        return
    
    def direct_force(self,mass,posn,soft):
        try:
            assert(self.nParticles > 0)
        except AssertionError:
            print("Please input data, no particles to compute!")
        try:
            assert(len(mass) == self.nParticles - 1)
        except AssertionError:
            print("Incorrect number of masses")
        try:
            assert(len(posn[0]) == self.nParticles - 1)
        except AssertionError:
            print("Incorrect number of positions")

        force = np.zeros((self.nParticles - 1,3))
        m2 = np.outer(mass,mass)
        if(np.isscalar(soft)):
            # Naive implementation, should be able to stop at nparticles/2?
            for i in range(self.nParticles-1):
                for j range(in self.nParticles-1):
                    if i == j:
                        continue
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (m2[i][j]) / (np.linalg.norm(posn[i] - posn[j]) + soft)**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
        else:
            for i in range(self.nParticles-1):
                for j range(in self.nParticles-1):
                    if i == j:
                        continue
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (m2[i][j]) / (np.linalg.norm(posn[i] - posn[j]) + soft[i])**3
                    force[i] += -1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
                    
        return force

class OctGrid(Gravity):
    def __init__(self,
                 filepath = None,
                 filename = "data.ascii"):
        if filepath is None:
            self.filepath = os.cwd()
        else:
            self.filepath = filepath
        self.filename = filename
        # NOTE Read in from file indexed from 1, not 0
        self.nParticles = -1
        self.nGasParticles = -1
        self.nStarParticles = -1

        self.particles = np.zeros(0)
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
        self.soft = s
    def set_pot(self,p):
        self.pot = p
        
    def set_variables():
        header,data = self.read_in_data(self.filepath,self.filename)
        self.nParticles = header[0]
        self.nGasParticles = header[1]
        self.nStarParticles = header[2]
        particles = np.zeros((self.nParticles))
        for i in range(self.nParticles):
            particles[i] = Particle([data[1][i],data[2][i],data[3][i]],
                                    [data[4][i],data[5][i],data[6][i]],
                                    data[0][i])
        self.soft = data[7]
        self.pot = data[8]
        self.particles = particles
        return

    def build_particles():
        particles = np.zeros((self.nParticles))
        for i in range(self.nParticles):
            particles[i] = Particle(self.posn[i],self.vels[i],self.mass[i])
        
    def close_enough(p1,p2):
        return

    def multipole(self, source, target):
        
        return

    def potential():
        return
    
    def compute_force():
        return
    
