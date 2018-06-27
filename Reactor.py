from CellList import CellList
from RuleBased import RuleBased
from Boundaries import Boundaries
from ParticleList import ParticleList
import Force
import numpy as np
from scipy.spatial.distance import cdist
from math import pi as PI
np.set_printoptions(threshold=np.nan)

#Particles in Format [amount, mass, radius, name]
# array with bindingSites[i] contains all binding sites of particleType i
# in spherical coordinates (phi, theta), referring to (azimuth, polar angle)
# max range: 0<phi<360, 0<theta<180
 
rawParticles    = [[100,35.0,5.0,"A"],[0,35.0,5.0,"B"],[0,35.0,5.0,"C"]] 
 
rawSites        = [{'a' : (PI, 0.0), 'b' : (PI, PI)},
                   {'d' : (PI, 0.0), 'e' : (PI, PI)},
                   {'g' : (PI, 0.0), 'h' : (PI, PI)}]

rawRules        = ["A(a) + A(b) -> A(a!1).B(b!1)"]
#-------------------------------------------------------------------------------

class Reactor:
    
    def __init__(self, volume, boundary=None, physics=None, forces=None):
        self.volume         = volume
        self.particles      = rawParticles
        self.n              = sum(particle[0] for particle in self.particles)
        self.b              = boundary

        self.particleList   = ParticleList(self.volume, self.particles, self.n, rawSites) # list with all particle objects

        self.maxRadius      = 2 * max(particle[2] for particle in self.particles)    # maximal radius of interactions
        self.cellList       = CellList(volume,self.maxRadius)                   # generates particle grid

        self.forceBond      = 0.5
        self.forceShell     = 0.5
        self.forceAngle     = 0.5
  
    def setup(self):
        self.pConfiguration = self.particleList.location_extracter()            # array containing all positions of the particles
        self.dConfiguration = self.particleList.diffusion_extracter()           # array containing all diffusions of i-th regarding particle        
        self.bConfiguration = self.cellList.add_particlelist(self.particleList) # array containing all box-indices of i-th particle
        self.nConfiguration = self.cellList.neighborlist_generator()            # array containing the list of neighbors of i-th Box
        self.boundaries     = Boundaries(self.n, self.volume)                   # Boundary conditions
        self.ruleNetwork    = RuleBased(self.particleList.pList, self.particles, rawRules)# empty interaction graph
        
    def apply_forces (self, stepsize, dT, disperse, componentDiffusion):
        l = self.n
        #Basic Force------------------------------------------------------------
        
        if disperse and dT < float("inf"):
            #print "BOOOM"
            allForce = Force.zero(l)
            compForce = np.sqrt(2 * dT * componentDiffusion) * Force.brownian(len(componentDiffusion))

            for c,component in enumerate(self.ruleNetwork.components):
                if len(component) > 1:
                    for i in component:
                        allForce[i] = compForce[c]
                else:
                    allForce[component[0]] = compForce[c]
        else:
            #g=5
            allForce = np.sqrt(2 * stepsize * self.dConfiguration) * Force.brownian(l)
        
        #allForce = Force.zero(l)
        #allForce = stepsize * Force.gravity(l)      
        
        #Additional Force-------------------------------------------------------
        allForce += Force.soft_shell(self, stepsize)
        #print Force.SoftShell(self,stepsize)
        allForce += Force.bonds(self, stepsize)
        #print Force.Bonds(self, stepsize)
        allForce += Force.angles(self, stepsize)

        #Adding all Forces------------------------------------------------------
        self.pConfiguration +=  allForce
        return 0

    def get_next_bireaction_time(self):
        # TODO currently only brute force O(#components^2) can be reduced to O(nlogn) or even O(n)!!!!
        # TODO components are also new created every timestep ... not updated
        pConfiguration     = self.pConfiguration
        
        components         = self.ruleNetwork.components
        componentDiffusion = self.ruleNetwork.componentDiffusion
        componentCenters   = []
        componentRadius    = []

        
        for component in components:
            if len(component) > 1:
                pConfigurationComponent   = pConfiguration[component]
                maxInnerComponentDistance = np.max(cdist(pConfigurationComponent,pConfigurationComponent))

                #would always work, but is quite slow for single particle component
                componentCenters.append(np.average(pConfigurationComponent,axis=0))
                componentRadius.append((maxInnerComponentDistance + self.maxRadius)/2.0)
            else:    
                componentCenters.append  (pConfiguration[component[0]])
                componentRadius.append   (self.particleList.pList[component[0]].radius)

        componentRadius    = np.array(componentRadius).reshape(len(componentRadius),1)       
        
        distances = cdist(componentCenters,componentCenters)
        np.fill_diagonal(distances, np.inf)
        
        def DistanceToHittime (distance, radius_c1, radius_c2, diff_c1, diff_c2):
            return ((distance - radius_c1 - radius_c2)/(np.sqrt(diff_c1) + np.sqrt(diff_c2)))**2/6.0
              
        times = DistanceToHittime(distances, componentRadius, componentRadius.transpose(), componentDiffusion, componentDiffusion.transpose())

        return np.min(times), componentDiffusion

    def location_infuser (self):
        pConfiguration = self.pConfiguration
        for i in range(len(pConfiguration)):
                self.particleList.pList[i].location = pConfiguration[i].tolist()
        return

    def __str__(self):
        r = ""
        for attr, value in self.__dict__.iteritems():
            if attr == "pList":
                r += "\n" + "ParticleList" + " : \n"
                for i in value:
                    r += str(i) + "\n"
            else:
                if "Config" not in `attr` and "cellList" not in `attr`:
                    r+= `attr` + " = " + str(value) + "\n"
        return r