import numpy as np

class Boundaries:
    def __init__(self, n, volume):
              
        self.minVol     = np.array([volume[::2]])
        self.maxVol     = np.array([volume[1::2]])
        self.tmp        = np.array([[-1,-1,-1]])
        
        self.tmp0       = np.repeat(self.minVol, n, axis = 0 )
        self.tmp1       = np.repeat(self.maxVol, n, axis = 0 )
        self.tmp2       = np.repeat(self.tmp,    n, axis = 0 )

    def reflective (self, reactor):
        pConfiguration  = reactor.pConfiguration
        
        ooLowerBound    = pConfiguration < self.tmp0                            # boolean array, if i-th dimension is exceeded
        ooUpperBound    = pConfiguration > self.tmp1                            # boolean array, if i-th dimension is unterschritten
        
        for component in reactor.ruleNetwork.components:
            if len(component) > 1:
                compConf            = pConfiguration[component]
                maxP                = np.max(compConf, axis=0)
                minP                = np.min(compConf, axis=0)
                
                leftLowerBoundary   = minP < self.minVol
                leftUpperBoundary   = maxP > self.maxVol        

                if leftLowerBoundary.any():
                    for particle in component:
                        pConfiguration[particle] -= (leftLowerBoundary * 2 * (minP - self.minVol)).reshape(3,)
                if leftUpperBoundary.any():
                    for particle in component:
                        pConfiguration[particle] -= (leftUpperBoundary * 2 * (maxP - self.maxVol)).reshape(3,)                
                
            else:
                particle = component[0]
                pLower = ooLowerBound[particle]
                pUpper = ooUpperBound[particle]
                
                if pLower.any():
                    pConfiguration[particle] -= pLower * 2 * (pConfiguration[particle] - self.minVol[0])
                if pUpper.any():
                    pConfiguration[particle] -= pUpper * 2 * (pConfiguration[particle] - self.maxVol[0])
                