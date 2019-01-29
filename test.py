import numpy as np
from gravity import *
from oct_tree import *

myTree = OctGrid("./","test_data.txt")
plist = myTree.buildParticleList()
myTree.buildTree(plist)
masses =[]
length = 0


for branch in myTree.tree.traverse():
    if branch.data:
        for d in branch.data:
            masses.append(d.pNum)

print(len(masses),masses)
