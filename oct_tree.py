import numpy as np
import sys
sys.setrecursionlimit(10000)

class Particle:
    def __init__(self,posn,velocity,mass,accel = None,error = None):
        self.posn = posn
        self.velocity = velocity
        self.mass = mass
        self.accel = None
        self.error = 0.0
        
        if accel is None:
            self.accel = np.zeros(3)
        else:
            self.accel = accel
        if error is not None:
            self.error = error
        return
    
    def addToAccel(self,a):
        self.accel += a.copy()
class Node:
    def __init__(self,posn,size,particles, parent = None):
        self.posn = posn
        self.size = size
        self.data = []
        self.nParticles = 0
        self.mass = 0.0
        if particles:
            self.data.append(particles)
            self.nParticles = 1
            self.mass = particles.mass
        self.CoM = None
        self.moments = None
        self.local_expansion = None
        self.parent = parent
        self.branches = [None,None,None,None,None,None,None,None]
        self.isLeaf = True
        return

    def setBranches(self,branches):
        self.branches = branches
    def __iter__(self):
        if not self.isLeaf:
            for branch in self.branches:
                yield branch
                
    def __len__(self):
        if self.data is not None:
            return len(self.data)


        
class OctTree():
    def __init__(self,size,maxNinNode = 1):
        self.root = self.addNode((0,0,0),size,None)
        self.size = size
        self.maxN = maxNinNode

        
    def buildTree(self,particles):
        for p in particles:
            self.insertParticle(self.root,self.size,self.root, p)
        return
    
    def addNode(self,posn,size,particle,parent = None):
        return Node(posn,size,particle,parent)

    def insertParticle(self, root, size, parent, particle):
        if root is None:
            new_posn = parent.posn
            shift = size/2.0
            branch = self.findBranch(parent,particle.posn)
            new_centre = (0.0, 0.0, 0.0)
            if branch == 0:
                # left down back
                new_centre = (new_posn[0] - shift, new_posn[1] - shift, new_posn[2] - shift )
                
            elif branch == 1:
                # left down forwards
                new_centre = (new_posn[0] - shift, new_posn[1] - shift, new_posn[2] + shift )
                
            elif branch == 2:
                # right down forwards
                new_centre = (new_posn[0] + shift, new_posn[1] - shift, new_posn[2] + shift )
                
            elif branch == 3:
                # right down back
                new_centre = (new_posn[0] + shift, new_posn[1] - shift, new_posn[2] - shift )

            elif branch == 4:
                # left up back
                new_centre = (new_posn[0] - shift, new_posn[1] + shift, new_posn[2] - shift )

            elif branch == 5:
                # left up forward
                new_centre = (new_posn[0] - shift, new_posn[1] + shift, new_posn[2] + shift )
                
            elif branch == 6:
                # right up forward
                new_centre = (new_posn[0] + shift, new_posn[1] - shift, new_posn[2] - shift )

            elif branch == 7:
                # right up back
                new_centre = (new_posn[0] + shift, new_posn[1] + shift, new_posn[2] - shift )
            return self.addNode(new_centre, size, particle, parent)
        elif not root.isLeaf:
            root.nParticles += 1
            root.mass += particle.mass
            branch = self.findBranch(root,particle.posn)
            new_size = root.size / 2.0
            root.branches[branch] = self.insertParticle(root.branches[branch], new_size, root, particle)
        elif root.isLeaf:
            if (len(root.data) + 1) <= self.maxN:
                root.mass += particle.mass
                root.data.append(particle)
                root.nParticles += 1
            else:
                root.data.append(particle)
                root.mass += particle.mass
                root.nParticles += 1
                dataList = root.data
                root.data = None
                root.isLeaf = False
                new_size = root.size/2.0
                #print(root.branches, dataList, dataList[0].posn, root.posn, root.size)
                for d in dataList:
                    branch = self.findBranch(root,d.posn)
                    root.branches[branch] = self.insertParticle(root.branches[branch], new_size, root, d)
        return root

    def checkPosition(self,node,posn):
        if posn[0] > node.posn[0] + node.size/2.0 or posn[0] < node.posn[0] - node.size/2.0:
            return False
        if posn[1] > node.posn[1] + node.size/2.0 or posn[1] < node.posn[1] - node.size/2.0:
            return False
        if posn[2] > node.posn[2] + node.size/2.0 or posn[2] < node.posn[2] - node.size/2.0:
            return False
        return True
    
    def getLeafAtPosn(self, root, posn):
        if root is None:
            return None
        elif root.isLeaf:
            return root.data
        else:
            branch = self.findBranch(root,posn)
            return self.findLeafAtPosn(root.branches[branch], posn)

    def findBranch(self, root, posn):
        # Magic Function for identifying which branch to pick
        DIRLOOKUP = {"3":0, "2":1, "-2":2, "-1":3, "1":4, "0":5, "-4":6, "-3":7}
        p1 = np.array(root.posn)
        p2 = np.array(posn)
        branch = 0
        for i in range(3):
            if p1[i] <= p2[i]:
                branch += (-4/ (i+1) / 2)
            else:
                branch += (4/ (i+1) / 2)
        branch = DIRLOOKUP[str(branch)]
        return branch

    def loopBranches(self,parent):
        i = 0
        for branch in parent.branches:
            if branch is None:
                continue
            if not branch.isLeaf:
                for subbranch in self.loopBranches(branch):
                    yield subbranch
            yield branch

    def traverse(self):
        if not self.root.isLeaf:
            for branch in self.loopBranches(self.root):
                yield branch


    
    
