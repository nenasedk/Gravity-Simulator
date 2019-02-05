from gravity import *
from oct_tree import *
from scipy.integrate import *

def dxdt(vel,dt):
    return vel

def dvdt(acc,dt):
    return acc

def EoM(w,t):
    posns = np.array([w[0],w[1],w[2]])
    vels = np.array([w[4],w[5],w[6]])
    acc = np.array([w[7],w[8],w[9]])
    dy2dt2 = dvdt

def IntegrateParticle(particle,tstep):
    particle.vels = odeint(dvdt,particle.accel,tstep)
    particle.posn = odeint(dvdt,particle.vels,tstep)
    return particle

def IntegrateTree(root,tstep):
    if root.isLeaf:
        for p in root.data:
            yield IntegrateParticle(p,tstep)
    else:
        for branch in root.branches:
            if branch is not None:
                if branch.isLeaf:
                    for p in branch.data:
                        yield IntegrateParticle(p,tstep)
                else:
                    IntegrateTree(branch,tstep)

def SolveTree(Octree,PartList,TStep,End,Order,Theta):
    time = np.linspace(0.0,End,TStep)
    for t in time:
        Octree.buildTree(Partlist)
        Octree.multipoleBuild(Octree.tree.root,Order)
        Octree.downwardExpansionTree(Octree.tree.root,Theta,Order)
        Octree.computeForce(Octree.tree.root,Theta,Order)
        Partlist = [p for p in IntegrateTree(Octree.tree.root,t)]
        #Plotting would go here
    return
                        
