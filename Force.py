import numpy as np
import sys
from math import sqrt, acos

#Basic Forces-------------------------------------------------------------------
def brownian (n):
    return np.random.normal(0,1,(n,3))

def gravity (n):
    return np.repeat([[9e10,-9e10,9e10]],n,axis=0)

def zero (n):
    return np.repeat([[0.0,0.0,0.0]],n,axis=0)
#-------------------------------------------------------------------------------  
  
  
#Custom Forces------------------------------------------------------------------
def soft_shell (reactor, stepsize):
    fs              = reactor.forceShell
    n               = reactor.n

    pList           = reactor.particleList.pList
    graph           = reactor.ruleNetwork.graph

    pConfiguration  = reactor.pConfiguration
    bConfiguration  = reactor.bConfiguration
    nConfiguration  = reactor.nConfiguration

    cellList        = reactor.cellList

    shellForce      = [[0.0,0.0,0.0] for _ in range(n)]

    for i in range(n):
        mass1, pos1, pre1, rad1 = pList[i].mass, pConfiguration[i], pList[i].forceprefactor, pList[i].radius       
        neighborbox = nConfiguration[bConfiguration[i]]
        
        for box in neighborbox:
            for neighborparticle in cellList.grid[box]: 
                j = neighborparticle
                if i != j:
                    force = 0.0
                    mass2, pos2, pre2, rad2 = pList[j].mass, pConfiguration[j], pList[j].forceprefactor, pList[j].radius
                    
                    dmax = rad1 + rad2
                    
                    direction = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]]
                    
                    dist = sqrt( direction[0]**2 + direction[1]**2 + direction[2]**2)
                    
                    if dist < dmax:
                        force += fs*(dmax - dist) / dist
                        
                        if pList[i].type in reactor.ruleNetwork.rules[pList[j].type]:
                            bond = reactor.ruleNetwork.check_for_interaction(pList[i], pList[j])
    
                            if bond and (i,j) not in graph.edges() and (j,i) not in graph.edges():
                                print "Bond formed between " + str(i) + " and " + str(j) + " at " + `bond`
                                graph.add_edge(i,j)
                                pList[i].sites[bond[0]][1] = j
                                pList[j].sites[bond[1]][1] = i
                                reactor.ruleNetwork.generate_angles()
                                reactor.ruleNetwork.set_componentproperties(reactor.particleList.pList)
           
                    shellForce[i][0] += stepsize * force * pre1 * direction[0] / mass1
                    shellForce[i][1] += stepsize * force * pre1 * direction[1] / mass1
                    shellForce[i][2] += stepsize * force * pre1 * direction[2] / mass1
                    
                    shellForce[j][0] -= stepsize * force * pre2 * direction[0] / mass2
                    shellForce[j][1] -= stepsize * force * pre2 * direction[1] / mass2
                    shellForce[j][2] -= stepsize * force * pre2 * direction[2] / mass2

    return shellForce


def bonds (reactor, stepsize):
    fb              = reactor.forceBond
    n               = reactor.n

    bonds           = reactor.ruleNetwork.graph.edges()
    pList           = reactor.particleList.pList
    pConfiguration  = reactor.pConfiguration

    bondForce       = [[0.0,0.0,0.0] for _ in range(n)]

    for (i,j) in bonds:
        force = 0.0
        mass1, pos1, pre1, rad1 = pList[i].mass, pConfiguration[i], pList[i].forceprefactor, pList[i].radius
        mass2, pos2, pre2, rad2 = pList[j].mass, pConfiguration[j], pList[j].forceprefactor, pList[j].radius
        
        dmax = rad1 + rad2
        
        direction = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]]        
        dist = sqrt( direction[0]**2 + direction[1]**2 + direction[2]**2)

        if dist > dmax:
            force = fb*(dist - dmax)/dist

        bondForce[i][0] -= stepsize * force * pre1 * direction[0] / mass1
        bondForce[i][1] -= stepsize * force * pre1 * direction[1] / mass1
        bondForce[i][2] -= stepsize * force * pre1 * direction[2] / mass1
        
        bondForce[j][0] += stepsize * force * pre2 * direction[0] / mass2
        bondForce[j][1] += stepsize * force * pre2 * direction[1] / mass2
        bondForce[j][2] += stepsize * force * pre2 * direction[2] / mass2
    return bondForce


def angles (reactor, stepsize):
    fa              = reactor.forceAngle
    n               = reactor.n
    
    angles          = reactor.ruleNetwork.angleList
    pList           = reactor.particleList.pList
    pConfiguration  = reactor.pConfiguration

    angleForce      = [[0.0,0.0,0.0] for _ in range(n)]

    for (i,j,k) in angles: # corresponds to the angle <JIK  -> 1-2-3
        
        mass1, pos1, pre1 = pList[j].mass, pConfiguration[j], pList[j].forceprefactor
        pos2              =                pConfiguration[i]
        mass3, pos3, pre3 = pList[k].mass, pConfiguration[k], pList[k].forceprefactor
        
        dir_j = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2]- pos2[2]]
        dir_k = [pos3[0] - pos2[0], pos3[1] - pos2[1], pos3[2]- pos2[2]]       
        
        dot_j  = dir_j[0] ** 2       + dir_j[1] ** 2       + dir_j[2] ** 2
        dot_jk = dir_j[0] * dir_k[0] + dir_j[1] * dir_k[1] + dir_j[2] * dir_k[2]  
        dot_k  = dir_k[0] ** 2       + dir_k[1] ** 2       + dir_k[2] ** 2

        # Calculate the direction of the force, middle particle will be shifted towards dir_j + dir_k
        # currently entire angle movement is based on the outside particles, middle is not considered
        # Two cases are important, 180 or any other angle
        
        nonZeroEntry = dir_k.index(filter(lambda x: x!=0, dir_k)[0])
        scalar = dir_j[nonZeroEntry]/dir_k[nonZeroEntry]
        
        #Case1: Directions are aligned, thus 180 degree
        if dir_k[0] * scalar == dir_j[0] and dir_k[1] * scalar == dir_j[1] and dir_k[2] * scalar == dir_j[2]  :
            #print "Case1"
            if dir_j[2] < dir_j[0]:
                mid = [dir_k[1],-dir_k[0],0.0]
            else:
                mid = [0.0, - dir_k[2], dir_k[1]]
            
        #Case2: Any other angle
        else:
            #print "Case3"
            mid = [dir_k[0] + dir_j[0], dir_k[1] + dir_j[1], dir_k[2] + dir_j[2] ]
        
        #Normalize the direction    
        lmid = sqrt(mid[0] ** 2 + mid[1]**2 + mid[2]**2)
        mid = [mid[0]/lmid, mid[1]/lmid, mid[2]/lmid]  
                        
        #Calculate the strength of the angle force
        len_j = sqrt(dot_j) 
        len_k = sqrt(dot_k)

        c = dot_jk/(len_j * len_k)
        if c > 1: c = 1
        if c < -1: c= -1
        
        angle = acos(c)
        
        # find angle0
        tmp = [pList[i].name]
        for site in pList[i].sites.keys():
            if pList[i].sites[site][1] in [j,k]:
                tmp += [site]
        if len(tmp) != 3:
            sys.exit("ERROR: multiple bonding ...")
        else:
            angle0 = reactor.ruleNetwork.rawAngles[tuple(tmp)]
            
        force = fa*(angle - angle0)**2
      
        if angle > angle0:
            force *= -1
        
        angleForce[j][0] -= stepsize * force * pre1 * mid[0] / mass1
        angleForce[j][1] -= stepsize * force * pre1 * mid[1] / mass1
        angleForce[j][2] -= stepsize * force * pre1 * mid[2] / mass1
        
        angleForce[k][0] -= stepsize * force * pre3 * mid[0] / mass3
        angleForce[k][1] -= stepsize * force * pre3 * mid[1] / mass3
        angleForce[k][2] -= stepsize * force * pre3 * mid[2] / mass3

        angleForce[i][0] += stepsize *force * mid[0] * (pre1 / mass1 + pre3 / mass3)
        angleForce[i][1] += stepsize *force * mid[1] * (pre1 / mass1 + pre3 / mass3)
        angleForce[i][2] += stepsize *force * mid[2] * (pre1 / mass1 + pre3 / mass3)

    return angleForce
