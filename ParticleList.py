from Particle import Particle
from RandomGenerator import RandomTripleUniform
from Tools import MassToDiffusion, MassToFactor, SitesToAngles
import numpy as np

class ParticleList:
    
    def __init__ (self, volume, particles, n, bindingSites):
        self.pList          = []
        self.n              = n
        pId                 = 0
        
        for i, particle in enumerate(particles):
            pType           = i + 1
            mass            = particle[1]
            radius          = particle[2]
            name            = particle[3]
            diffusion       = MassToDiffusion(radius)
            sites           = {site: [bindingSites[i][site],-1] for site in bindingSites[i]}
            angles          = SitesToAngles(sites)
    
            forceprefactor  = MassToFactor(radius)
    
            #TESTCASES---------------------------
            #seeds = [[[0.0,130.0,4.0], [30.0,125.0,4.0],[60.0,115.3,4.0], [90.0, 93.0, 4.0],[110.0,70.0,4.0],[120.0,50.0,4.0],[130.0,0.0,4.0]], [[2.0,0.0,4.0]]]        
            #seeds = [[[0.0,0.0,0.0], [0.0, 9.0,00.0], [0.0, 18.0,00.0]]]
            #seeds = [[[0.0,0.0,0.0],[9.0,0.0,0.0]]]
            
            for _ in range(particle[0]):
                seed        = RandomTripleUniform(volume)
                #seed        = seeds[i][_]
                self.pList += [Particle(pId, pType, seed, mass, radius, name, diffusion, sites, angles, forceprefactor)]
                pId        += 1
                
    def location_extracter (self):
        return np.array( [self.pList[i].location   for i in range(self.n)])

    def diffusion_extracter (self):
        return np.array([[self.pList[i].diffusion] for i in range(self.n)])
