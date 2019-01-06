import numpy as np
class Particle:
    def __init__(self,posn,velocity,mass):
        self.posn = posn
        self.velocity = velocity
        self.mass = mass
        return
    
class Node:
    def __init__(self,posn,size,data):
        self.posn = posn
        self.size = size
        self.data = data

        self.branches = [None,None,None,None,None,None,None,None]
        self.isLeaf = True
        return
    
    
class OctTree():
    def __init__(self,size,maxNinNode = 10):
        self.root = self.add_node((0,0,0),size,[])
        self.size = size
        self.maxN = maxNinNode

    def addNode(self,posn,size,data):
        return Node(posn,size,data)

    def insertNode(self, root, size, parent, data):
        if root is None:
            new_posn = parent.posn
            shift = size/2.0
            branch = self.findBranch(parent,data.posn)
            new_centre = (0.0, 0.0, 0.0)
                        if branch == 0:
                # left down back
                new_centre = (new_posn[0] - offset, new_posn[1] - offset, new_posn[2] - offset )
                
            elif branch == 1:
                # left down forwards
                new_centre = (new_posn[0] - offset, new_posn[1] - offset, new_posn[2] + offset )
                
            elif branch == 2:
                # right down forwards
                new_centre = (new_posn[0] + offset, new_posn[1] - offset, new_posn[2] + offset )
                
            elif branch == 3:
                # right down back
                new_centre = (new_posn[0] + offset, new_posn[1] - offset, new_posn[2] - offset )

            elif branch == 4:
                # left up back
                new_centre = (new_posn[0] - offset, new_posn[1] + offset, new_posn[2] - offset )

            elif branch == 5:
                # left up forward
                new_centre = (new_posn[0] - offset, new_posn[1] + offset, new_posn[2] + offset )
                
            elif branch == 6:
                # right up forward
                new_centre = (new_posn[0] + offset, new_posn[1] - offset, new_posn[2] - offset )

            elif branch == 7:
                # right up back
                new_centre = (new_posn[0] + offset, new_posn[1] + offset, new_posn[2] - offset )
            return self.addNode(new_centre, size, data)
        elif root.posn != data.posn:
            if root.isLeaf:
                if (len(root.data) + len(data.data)) < self.maxN:
                    root.data.append(data)
                else:
                    root.data.append(data)
                    root.isLeaf = False
                    dataList = root.data
                    new_size = root.size/2.0
                    for d in data:
                        branch = self.findBranch(root,d.posn)
                        root.branches[branch] = self.insertNode(root.branches[branch], new_size, root, d)
            else:
                branch = self.findBranch(root,data.posn)
                new_size = root.size / 2.0
                root.branches[branch] = self.insertNode(root.branches[branch], new_size, root, data)
        return root


    def getLeafAtPosn(self, root, posn):
        if root is None:
            return None
        elif root.isLeaf:
            return root.data
        else:
            branch = self.findBranch(root,posn)
            return self.findLeafAtPosn(root.branches[branch], posn)

    def findBranch(self, root, posn):
        p1 = root.posn
        p2 = posn
        d = p2 - p1
        branch = 0
        for i,x in enumerate(d):
            if x <= 0:
                branch += (-4/ (i+1) / 2)
            else:
                branch += (4/ (i+1) / 2)
        return branch

    
