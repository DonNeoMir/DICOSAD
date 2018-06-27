import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations, product
from Tools import MassToDiffusion
import numpy as np
import sys
from math import acos, cos, sin

class RuleBased(object):

    def __init__(self, pList, particles, rawRules):
        n = len(pList)
        self.nodes      = range(n)
        self.graph      = nx.Graph()
        self.graph.add_nodes_from(self.nodes)
        self.components = []
        self.set_componentproperties(pList)
        
        self.rawAngles  = self.create_angles(pList)
        self.angleList  = []
        print self.rawAngles
        self.rawRules   = rawRules
        self.rules      = self.create_rules(particles)
        
        #self.graph.add_edge(0,1)
        #self.graph.add_edge(1,2)
        #self.graph.add_edge(2,3)

    def visualize(self):
        nx.draw(self.Graph)
        plt.show()

    def create_angles(self,pList):
        rawAngles = {}
        
        for particle in pList:
            for (i,j) in list(combinations(particle.sites.keys(),2)):
                if (particle.name,i,j) not in rawAngles:
                    att1 = particle.sites[i][0][0]
                    lon1 = particle.sites[i][0][1]
                    
                    att2 = particle.sites[j][0][0]
                    lon2 = particle.sites[j][0][1]
                    
                    angle = acos(sin(att1)*sin(att2) + cos(att1)*cos(att2)*cos(lon1-lon2))
                    rawAngles[(particle.name,i,j)] = angle
                    rawAngles[(particle.name,j,i)] = angle        
        return rawAngles

    def generate_angles(self):
        for node in self.nodes:
            neighbors = list(nx.all_neighbors(self.graph,node))
            angles = [[node] + list(pair) for pair in (combinations(neighbors,2))]
            for angle in angles:
                if angle not in self.angleList:
                    self.angleList += [angle]
                    
                        
    def create_rules(self,particles):
        rules = {}
        nameToType = {}
        for i,particle in enumerate(particles):
            name = particle[3]
            pType = i + 1
            nameToType[name] = pType
            if pType in rules:
                sys.exit("ERROR: Double particle name!")
            else:
                rules[pType] = {}
                
        for rule in self.rawRules:
            left = rule.split("->")[0].split()
                        
            edu1, site1 = left[0].split("(")
            edu2, site2 = left[2].split("(")
            
            site1 = site1[:-1]
            site2 = site2[:-1]
            
            if edu1 not in nameToType or edu2 not in nameToType:
                sys.exit("ERROR: Invalid educt!")
            else:
                edu1 = nameToType[edu1]
                edu2 = nameToType[edu2]
            if edu2 in rules[edu1]:
                if (site1,site2) not in rules[edu1][edu2]:
                    rules[edu1][edu2] += [(site1,site2)]
            else:
                rules[edu1][edu2] = [(site1,site2)]
                
            if edu1 in rules[edu2]:
                if (site2,site1) not in rules[edu2][edu1]:
                    rules[edu2][edu1] += [(site2,site1)]
            else:
                rules[edu2][edu1] = [(site2,site1)] 
        return rules
    
    def check_for_interaction(self, particle_i, particle_j):
        rules = self.rules
        freeSitesI = [i for i in particle_i.sites if particle_i.sites[i][1] == -1]
        freeSitesJ = [i for i in particle_j.sites if particle_j.sites[i][1] == -1]
        freeInteractions = list(product(freeSitesI,freeSitesJ))

        for inter in freeInteractions:
            if inter in rules[particle_i.type][particle_j.type]:
                return inter
        return 0
        
        
    def set_componentproperties(self, pList):
        self.components = list(nx.connected_components(self.graph))
        
        componentDiffusion = []        
        
        for component in self.components:
            if len(component) > 1:
                massComponent = sum(pList[i].mass for i in component)

                componentDiffusion.append(MassToDiffusion(massComponent))
            else:    
                componentDiffusion.append(pList[component[0]].diffusion)

        self.componentDiffusion = np.array(componentDiffusion).reshape(len(componentDiffusion),1)    
            
    def __str__(self):
        return `self.Graph.edges()`
    
        