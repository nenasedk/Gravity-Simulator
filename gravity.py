from oct_tree import *
import numpy as np
from numpy.linalg import norm
from scipy import special
import os 
import sys
sys.setrecursionlimit(2000)

class Gravity:
    def __init__(self):
        return

    def read_in_data(self,filepath, filename):
        if not filepath.endswith('/'):
            filepath += '/'
        if not filename.endswith('.txt'):
            if not filename.endswith('.ascii'):
                print("Can only read in .txt or .ascii data files.")
        
        with open(filepath + filename) as f:
            header = list(map(int,f.readline().split()))
            f.close()
        data = np.genfromtxt(filepath + filename,skip_header=1)
        data = data.reshape(9,header[0])
        return header, data

    def density(self,mass,posn,n,log = False,consth = False):
        coords = self.cart_to_sphere(posn)
        coordinds = coords[:,0].argsort()
        s_coords = coords[:,0][coordinds]

        s_mass = mass[coordinds]
        max_rad = np.max(s_coords)

        lbase = np.e
        if consth:
            length = int(np.ceil(self.nParticles/n))
            radout = np.zeros(length)
            density = np.zeros(length)
            error = np.zeros(length)
            radius = 0.0
            j = 0
            for i in range(length):
                start = j
                stop = int(np.min([j+n,self.nParticles-1]))
                volume = (4.0/3.0)*np.pi * (s_coords[stop] - s_coords[start])**3.0
                m_shell = 0.0
                while(j<=stop and (j+1)<self.nParticles):
                    m_shell += s_mass[j]
                    j+=1
                density[i] = m_shell/volume
                radout[i] = (s_coords[stop] + s_coords[start])/2.0
                error[i] = np.sqrt(j - start)
            return radout,density,error,s_coords
        radout = np.zeros(n-1)
        density = np.zeros(n-1)
        error = np.zeros(n-1)
        radius = 0.0
        if log:
            h = np.logspace(np.log(np.min(s_coords)),
                            np.log(max_rad),
                            n,base = lbase)
            radius = 0.0
            j = 0
            for i in range(n-1):
                volume = (4.0/3.0)*np.pi * (h[i+1]- h[i])**3
                m_shell = 0.0
                errmass = []
                while(s_coords[j]<h[i+1] and j<self.nParticles-1):
                    m_shell += s_mass[j]
                    errmass.append(s_mass[j])
                    j+=1
                density[i] = m_shell/volume
                radout[i] = np.power(10,(np.log10(h[i+1]) + np.log10(h[i]))/2.0)
                error[i] = np.sqrt(len(errmass))
            return radout,density,error,s_coords 
        else:
            h = max_rad/n
            radius = 0.0
            j = 0
            for i in range(n-1):
                volume = (4.0/3.0)*np.pi * (((i+1)*h) - (i*h))**3
                m_shell = 0.0
                errmass = []
                while(s_coords[j]<(i+1)*h):
                    m_shell += s_mass[j]
                    errmass.append(s_mass[j])
                    j+=1
                density[i] = m_shell/volume
                radout[i] = (i*h + (i+1)*h)/2.0
                error[i] = np.sqrt(len(errmass))
            return radout,density,error,s_coords         
            
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
        s_coords = coords[:,0][coordinds]
        s_mass = mass[coordinds]
        m = 0.0
        i = 0 
        while(m<M_t/2): # Find Half-Mass Radius
            m += s_mass[i]
            i += 1            
        r_half = s_coords[i]
        V_p = np.sqrt(M_t/(2.0 * r_half)) # Velocity
        return r_half/V_p # Return crossing time.
    
    def cart_to_sphere(self,data):
        if len(data) < 4:
            dist = np.square(data[0]) + np.square(data[1]) + np.square(data[2])
            r = np.sqrt(np.square(data[0]) + np.square(data[1]) + np.square(data[2]))
            tht = np.arctan2(data[1],data[0])
            phi = np.arccos(data[2]/r)
            newcoords = (r,tht,phi)
        else:
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
            assert(len(mass) == len(posn))
        except AssertionError:
            print("Incorrect number of masses or positions")
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
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j])**2.0 + soft**2.0)**(1.5)
                    force[i] = -1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
        else:
            for i in range(self.nParticles-1):
                for j in range(self.nParticles-1):
                    if i == j:
                        continue
 
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j])**2/0 + soft[i]**2.0)**(1.5)
                    force[i] += 1.0*(posn[i] - posn[j]) * scalar
                    #force[j] = force[j] - force[i]
                    
        return force

    """
    Uses Newton's third law to reduce number of computations
    """
    def direct_force_fast(self,mass,posn,soft):
        try:
            assert(self.nParticles > 0)
        except AssertionError:
            print("Please input data, no particles to compute!")
            sys.exit(1)
        try:
            assert(len(mass) == len(posn))
        except AssertionError:
            print("Incorrect number of masses")
            sys.exit(1)

        force = np.zeros((self.nParticles ,3))
        # m2 = np.outer(mass,mass)
        if(np.isscalar(soft)):
            # Naive implementation, should be able to stop at nparticles/2?
            for i in range(0,self.nParticles-1):
                for j in range(self.nParticles-1,i,-1):
                    if i == j:
                        break
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j])**2 + soft**2.0)**(1.5)
                    fij = (posn[i] - posn[j]) * scalar
                    force[i] += fij
                    force[j] += -1.0 * fij
        else:
            for i in range(0,self.nParticles-1):
                if(i%500 == 0):
                    print(str(i/500) + "%")
                for j in range(self.nParticles-1,i,-1):
                    if i == j:
                        break
                    # Should check if softening is a single scalar value, or unique for each element
                    scalar = (mass[i]*mass[j]) / (np.linalg.norm(posn[i] - posn[j])**2.0 + soft[i]**2.0)**(1.5)
                    fij = (posn[i] - posn[j]) * scalar
                    force[i] += fij
                    force[j] += -1.0 * fij
        return force









"""
Implementation of a fast multipole method based off of Greenwood et al.

"""    
class FMM(Gravity):
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
        
    def buildParticleList(self,pNum = None):
        if pNum:
            header,data = self.read_in_data(self.filepath,self.filename)
            self.nParticles = pNum
            self.nGasParticles = int(header[1])
            self.nStarParticles = pNum
            particles = []
            self.size = np.max([np.max(data[1][0:pNum]),
                                np.max(data[2][0:pNum]),
                                np.max(data[3][0:pNum])])*2.0
            for i in range(self.nParticles):
                particles.append(Particle(np.array([data[1][i],data[2][i],data[3][i]]),
                                          np.array([data[4][i],data[5][i],data[6][i]]),
                                          data[0][i]))
            self.soft = data[7][0:pNum]
            self.s_soft = np.mean(data[7][0:pNum]) *10
            self.pot = data[8][0:pNum]
        else:
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
            self.s_soft = np.mean(data[7])
            self.pot = data[8]
        return particles

    """
    Tree Structure
    """
    def buildTree(self,particles):
        self.tree = OctTree(self.size,1)
        self.tree.buildTree(particles)
        return
    
    def distanceCriterion(self,posn,node,eps):
        if node.isLeaf:
            return False
        if self.tree.checkPosition(node,posn):
            return False
        #print (node.size/(np.linalg.norm(particle.posn - node.posn)))
        if(node.size/(np.linalg.norm(np.array(posn) - np.array(node.posn))) < eps):
            return True
        else:
            return False

    """ 
    Multipole Functions
    Upward Steps
    """
    def multipoleBuild(self,node,L):
        if node.isLeaf:
            self.calcMultipoleMoments(node,L)
            return node
        else:
            node.moments = {}
            for branch in node.branches:
                if branch is None:
                    continue
                #self.multipoleBuild(branch,node,L,M)
                cnode = self.multipoleBuild(branch,L)
                #moments.append(self.shiftMultipoleMomement(branch.posn,node.posn, branch.moments))
                self.shiftMultipoleMoments(node,cnode,L)
                # multipole expansion to local expansion
            return node
            
    def calcMultipoleMoments(self,node,L):
        node.moments = {}
        for l in range(0,L+1):
            for m in range(-l,l+1):
                node.moments[(l,m)] = 0.0
        for l in range(L+1):
            phi_l = 0
            for m in range(-l,l+1):
                node.moments[(l,m)] = self.calcMultipoleCoefficient(node,l,m)        
        return node
    
    def calcMultipoleCoefficient(self,node,l,m):
        Q = 0
        for p in node.data:
            r2 = self.cart_to_sphere(p.posn)
            sph = special.sph_harm(m,l,r2[1],r2[2])
            Q += p.mass * r2[0]**l * np.conj(sph)
        return Q.real
    
    def shiftMultipoleMoments(self, parent, child, L):
        #shift_moments = {}
        if parent.moments == {}:
            for l in range(0,L+1):
                for m in range(-l,l+1):
                    parent.moments[(l,m)] = 0.0
                    
        polar = self.cart_to_sphere(np.array(child.posn))
        for j in range(0,L+1):
            for k in range(-j,j+1):
                M_lm = 0.0
                for n in range(0,j+1):
                    for m in range(-n,n+1):
                        try:
                            M_lm += (child.moments[(j-n,k-m)] * 1j**(np.abs(k) - np.abs(m) - np.abs(k-m)) \
                                     * self.expansionCoeff(n,m)*self.expansionCoeff(j-n,k-m) \
                                     * polar[0]**n * special.sph_harm(-m,n,polar[1],polar[2]))\
                                     /self.expansionCoeff(j,k)
                        except:
                            continue
                parent.moments[(j,k)] = M_lm
                #print parent.moments[(j,k)]
        return

    def expansionCoeff(self,n,m):
        return ((-1)**n) / np.sqrt(np.math.factorial(n-m)*np.math.factorial(n+m))

    """
    Local Expansions
    Downward Steps
    """        
    def multipoleExpand(self,target,nodelist,eps,L):
        nlist = []

        if target.isLeaf and nodelist == []:
            #print target.data[0].accel,target.data[0].posn
            return
        
        for node in nodelist:
            if node is None:
                continue
            elif self.posnComp(node.posn,target.posn):
                for b in target.branches:
                    if b is not None:
                        nlist.append(b)
                continue

            elif target.isLeaf and node.isLeaf:
                # This is super ugly because I was having issues
                # with all accels actually being the same array
                target.data[0].addToAccel(self.directForce(node.data[0],target.data[0]))
                continue
            
            elif self.distanceCriterion(target.posn,node,eps):
                self.shiftLocalExpansion(target,node,L)
                continue

            else:
                for b in node.branches:
                    if b is not None:
                        nlist.append(b)

        if target.isLeaf:
            self.multipoleExpand(target,nlist,eps,L)
                        
        else:
            for tbranch in target.branches:
                if tbranch is not None:
                    self.shiftLocalExpansion(tbranch,target,L)
                    self.multipoleExpand(tbranch,nlist,eps,L)

                
    def localExpansionMoment(self,target,node, L):
        polar = self.cart_to_sphere(node.posn)
        if not target.local_expansion:
            target.local_expansion = {}
            for l in range(0,L+1):
                for m in range(-l,l+1):
                    target.local_expansion[(l,m)] = 0.0
                    
        for j in range(0,L+1):
            for k in range(-j,j+1):
                L_kj = 0.0
                for n in range(0,L+1):
                    for m in range(-n,n):
                        try:
                            L_kj += (node.moments[(n,m)]*1j**((np.abs(k-m)-np.abs(k) - np.abs(m))) \
                                    *self.expansionCoeff(n,m) * self.expansionCoeff(j,k) \
                                    * special.sph_harm(m-k,j+n,polar[1],polar[2]))\
                                    /((-1)**n * self.expansionCoeff(j+n,m-k) * polar[0]**(j+n+1))
                        except:
                            continue
                target.local_expansion[(j,k)] += L_kj
        return 

    def shiftLocalExpansion(self,target,node,L):
        #polar = self.cart_to_sphere(np.array(target.posn) - np.array(node.posn))
        polar = self.cart_to_sphere(np.array(node.posn))
        if not target.local_expansion:
            target.local_expansion = {}
            for l in range(0,L+1):
                for m in range(-l,l+1):
                    target.local_expansion[(l,m)] = 0.0      
        for j in range(0,L+1):
            for k in range(-j,j+1):
                L_kj = 0.0
                for n in range(j,L+1):
                    for m in range(-n,n+1):
                        try:
                            top = node.moments[(n,m)]*1.0j**((np.abs(m)-np.abs(m-k) - np.abs(k))) \
                                  * self.expansionCoeff(n-j,m-k) * self.expansionCoeff(j,k) \
                                  * special.sph_harm(m-k,n-j,polar[1],polar[2])*polar[0]**(n-j)
                            bot = ((-1)**(n+j) * self.expansionCoeff(n,m))
                            L_kj += top/bot
                        except:
                            continue
                target.local_expansion[(j,k)] += L_kj
        return 

    """
    Potential Evaluation
    """
    def potentialExpansion(self,node,evalPosn, L):
        p2 = self.cart_to_sphere(evalPosn)
        pot = {}
        for l in range(0,L+1):
            for m in range(-l,l+1):
                pot[(l,m)] = node.local_expansion[(l,m)] * special.sph_harm(m,l,p2[1],p2[2]) * p2[0]**l
        return pot
    
    def computeForce(self,node,eps,L):
        if node.isLeaf:
            for p in node.data:
                #print p.accel
                pot = self.potentialExpansion(node,p.posn,L)
                #print pot[(0,0)]
                p.addToAccel(-1.0* np.array([np.real(pot[(1,1)]),np.imag(pot[(1,1)]),np.real(pot[(1,0)])]))
            return
        else:
            for branch in node.branches:
                if branch is not None:
                    self.computeForce(branch,eps,L)

    def directForce(self,p1,p2):
        fij = np.zeros(3)
        if self.posnComp(p1.posn,p2.posn):
            return fij
        scalar = p1.mass / (np.linalg.norm(np.array(p2.posn).copy() - np.array(p1.posn).copy())**2.0 + self.s_soft**2.0)**(1.5)
        fij = -1.0*(np.array(p2.posn).copy() - np.array(p1.posn).copy()) * scalar
        return fij.copy()
    
    def downwardExpansionTree(self,root,eps,L):
        for branch in root.branches:
            if branch is not None:
                self.multipoleExpand(branch,root.branches,eps,L)
                

                
    def posnComp(self,tuple1,tuple2):
        if tuple1[0] == tuple2[0] and tuple1[1] == tuple2[1] and tuple1[2] == tuple2[2]:
            return True
        else:
            return False
