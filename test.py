import numpy as np
from gravity import *
from oct_tree import *
import matplotlib.pyplot as plt
ORDER = 4
THETA = 0.01
myTree = FMM("./","data.ascii")
plist = myTree.buildParticleList(5000)
myTree.buildTree(plist)
myTree.multipoleBuild(myTree.tree.root,ORDER)
myTree.downwardExpansionTree(myTree.tree.root,THETA,ORDER)
myTree.computeForce(myTree.tree.root,THETA,ORDER)
forces = []
masses = []
posns = []
print "\n\n\n\n"
def looptree(root):
    for branch in root.branches:
        #print branch.local_expansion
        if branch is None:
            continue
        if branch.isLeaf:
            #print branch.local_expansion
            #for m in branch.local_expansion:
            #    pot.append(branch.local_expansion[m])      
            for p in branch.data:
                #print p.accel
                masses.append(p.mass)
                posns.append(p.posn)
                forces.append(p.accel)
        else:
            if branch is not None:
                looptree(branch)
looptree(myTree.tree.root)
fig,ax = plt.subplots(figsize = (12,8))
polar = myTree.cart_to_sphere(np.array(posns))

print np.mean(np.array(polar[:,0]))
#print posns np.linalg.norm(np.array(forces),axis = 1)*np.array(masses)
ax.scatter(polar[:,0], np.linalg.norm(np.array(forces),axis = 1)*np.array(masses),marker = '.',s = 15,linewidth = 2, c = u'#ff7f0e')
ax.set_title("Radial Force distribution using order " + str(ORDER) + " FMM, with $ \Theta $ = " + str(THETA))
ax.set_xlabel("Radius")
ax.set_ylabel("Force")
ax.set_xlim(0,20)
ax.set_ylim(0,2e6)
#ax.set_yscale("log")
plt.show()
