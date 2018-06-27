import copy

class Particle:

    def __init__(self, idP, typeP, initial, mass, radius, name, diffusion, sites, angles, forceprefactor):
        
        self.id             = idP
        self.type           = typeP
        self.name           = name
        
        self.location       = initial
        self.initial        = initial
        
        self.mass           = mass
        self.radius         = radius
        self.diffusion      = diffusion
        self.sites          = copy.deepcopy(sites)
        self.angles         = angles

        self.forceprefactor = forceprefactor
        #Temporary for bindings ..........
        self.bondList       = 0

    def __str__(self):
        r = ''
        for attr, value in self.__dict__.iteritems():
            r+= `attr` + ' = ' + str(value) + "\n"
        return r
    
    def get_LAMMPS(self):
        r = str(self.id) + " " + str(self.type) + " " + str(self.location[0]) + " " + str(self.location[1]) + " " + str(self.location[2])
        return r