import numpy as np
from gravity import *
from oct_tree import *

myTree = FMM("./","test_data.txt")
plist = myTree.buildParticleList()
myTree.buildTree(plist)
myTree.multipoleBuild(myTree.tree.root,2,2)
myTree.multipoleExpand(myTree.tree.root,2,2)
myTree.computeForce(myTree.tree.root,2,2)
masses =[]
length = 0


for branch in myTree.tree.traverse():
    if branch.isLeaf:
        for p in branch.data:
            masses.append(p.accel)

print(np.array(masses))
#print(np.array(masses)[:,1])
#print(np.array(masses)[:,2])

