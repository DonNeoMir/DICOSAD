from timeit import default_timer as timer
#time measurement not that cool yet, measures wall clock not cpu....

from Simulation import Simulation

def main(timesteps = 200000):
    #print 'Start!'
    start = timer()

    s = Simulation(timesteps)
    s.setup()
    time = s.start_simulation()

    #s.reactor.RuleNetwork.Visualize()
    c = s.reactor.ruleNetwork.components

    #print 'Done!'
    end = timer()
    #print "CPU seconds elapsed:    " + str(end - start) 
    #print "System seconds elapsed: " + str(time)  
    #print sum(map(lambda x: len(x),c))/(len(c)+0.0)
    #print c
    return c
    #return [sum(map(lambda x: len(x),c))/(len(c)+0.0), end - start, time]

if __name__ == '__main__':
    main()
