import numpy as np
import math
from itertools import product as Prod
from matplotlib.pyplot import box

class CellList:
    
    def __init__(self, volume, r_cut):
        dim = np.array([volume[i+1] - volume[i] for i in range(0,len(volume),2)])

        self.r_cut = r_cut

        self.origin = volume[::2]
        self.maxd = self.origin + dim
        self.blocks = map(int,np.ceil(dim/r_cut))
       
        xn, yn, zn = self.blocks
        self.grid = [[] for _ in range (xn * yn * zn) ]

    def __str__(self):
        return `self.grid`


    def add_particle(self,particle):
        box = self.get_box(particle, particle.location)

        self.grid[box] += [particle.id]
        
        return box

    def add_particlelist (self, particleList):
        return [self.add_particle(particleList.pList[i]) for i in range(particleList.n)]
    
    
    def update_particle(self, particle, box, pos):
        idP = particle.id

        newBox = self.get_box(particle, pos)
        
        if newBox != box:
            if newBox < 0 or newBox >= len(self.grid):
                print "ERROR: particle out of bounds"
                print idP, pos, box, newBox, len(self.grid)
                exit()
            self.grid[box].remove(idP)
            self.grid[newBox] += [idP]
            return newBox

        return -1
    

    def update_particlelist(self,reactor):
        pList = reactor.particleList.pList
        bConfiguration = reactor.bConfiguration
        pConfiguration = reactor.pConfiguration
        for i in range(len(pList)):
            box = self.update_particle(pList[i], bConfiguration[i], pConfiguration[i])
            if box > -1:
                bConfiguration[i] = box
        
    def get_box(self, particle, pos):
        
        r_cut = self.r_cut
        origin = self.origin
        bx, by, bz = self.blocks
        
        ix = int(math.floor((pos[0] - origin[0])/r_cut))
        iy = int(math.floor((pos[1] - origin[1])/r_cut))
        iz = int(math.floor((pos[2] - origin[2])/r_cut))
        
        if ix == bx : ix -=1
        if iy == by : iy -=1
        if iz == bz : iz -=1
        return  (ix*by + iy)*bz + iz
        
 
    def neighborlist_generator(self):
        xn, yn, zn = self.blocks
        neighborIds = [[] for _ in range (xn * yn * zn) ]
        
        allNeighbors = list(Prod([-1,0,1],repeat = 3))

        for x in range(xn):
            for y in range(yn):
                for z in range(zn):
                    middlebox = int((x*yn + y)*zn + z)

                    for seed in allNeighbors:
                        ix = x + seed[0]
                        iy = y + seed[1]
                        iz = z + seed[2]
                        
                        box = (ix*yn + iy)*zn + iz

                        if -1 < ix < xn  and -1 < iy < yn  and -1 < iz < zn :
                            if middlebox in neighborIds[box]:
                                continue
                            neighborIds[middlebox] += [box]
        return neighborIds