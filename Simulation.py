from Reactor import Reactor
import Output

class Simulation:  
    
    def __init__(self, timesteps):
        #variables should been passed from files ... but not yet!!
        self.name      = 'Trials1'
        
        self.timesteps = timesteps
        self.stepsize  = 1e-9
        self.threshold = 150#float("inf")

        self.volume    = [-300,300,-300,300,-300,300]      
        
        boundary       = None
        physics        = None
        forces         = None

        self.reactor   = Reactor(self.volume, boundary, physics, forces)  
    
    def setup(self):
        Output.NewLAMMPSFile(self.name)
        Output.NewMSDFile   (self.name)
        Output.NewTCLScript (self.name,self.reactor)
        self.reactor.setup()
        return
    
    def start_simulation(self):
        time = 0.0
        Output.DumpOneTimeStep(self.name,self.reactor,0)
        
        while time < 0.05:
        #for i in range(1,self.timesteps + 1):

            dT, componentDiffusion = self.reactor.get_next_bireaction_time()
            disperse = dT/self.stepsize > self.threshold
            self.reactor.apply_forces(self.stepsize, dT, disperse, componentDiffusion)
            
            #if not i % 100:
            #    print i
            #    self.reactor.location_infuser()
            #    Output.DumpOneTimeStep(self.name,self.reactor,i)

            self.reactor.boundaries.reflective(self.reactor)         
            self.reactor.cellList.update_particlelist(self.reactor)

            time += dT if disperse else self.stepsize
            
        #print time
        Output.NewMSDScript(self.name)
        return time