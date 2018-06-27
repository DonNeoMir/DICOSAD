import math as ma

Boltzmann = 1.38064852e-23 #J K^-1

def MassToDiffusion (R):
    # gets a mass in kDa and converts it to its corresponding Diffusion Coefficient
    # gets a radius in nm and converts it to its corresponding Diffusion Coefficient
    
    #if only mass is given, use formula below------------------
    #R = 0.066*(m*1e3)**(1/3.) #mass kDa -> nm radius of sphere
    #----------------------------------------------------------
    
    D = (300*Boltzmann)/(6*ma.pi*6.75*R*1e-12) # radius of sphere nm-> Diffusion m^2/s    
    #output in nm^2/s = 1e-6 um^2/s
    return D*1e18

def MassToFactor (R):
    #calculates the weird prefactor of particle ... compare readdy paper
    return 1/(6*ma.pi*6.75*R*1e-12)


def SitesToAngles(sites):
    from itertools import combinations
    angles      = {}
    SiteNames   = [a for a in sites]

    for pair in combinations(SiteNames,2):
        phi1    = sites[pair[0]][0][0]
        phi2    = sites[pair[1]][0][0]
        
        theta1  = sites[pair[0]][0][1]
        theta2  = sites[pair[1]][0][1]
        
        angle   = ma.acos( ma.sin(theta1 + ma.pi) * ma.sin(theta2 + ma.pi) + ma.cos(theta1 + ma.pi) * ma.cos(theta2 + ma.pi) * ma.cos(abs(phi1 - phi2)))
        
        angles[tuple(sorted(pair))] = angle
    return angles


