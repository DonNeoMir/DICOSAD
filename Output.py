def NewLAMMPSFile(name):
    open('LAMMPS_{}.lammpstrj'.format(name),'w+')
 
def DumpOneTimeStep(name,reactor, timestep):
    pList    = reactor.particleList.pList
    n        = reactor.n
    vol      = map(str,reactor.volume)
    
    f = open('LAMMPS_{}.lammpstrj'.format(name),'a')
    f.write("ITEM: TIMESTEP\n")
    f.write(`timestep` + "\n")
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write(`n` + "\n")
    f.write("ITEM: BOX BOUNDS\n")
    f.write(vol[0] + " " + vol[1] + "\n")
    f.write(vol[2] + " " + vol[3] + "\n")
    f.write(vol[4] + " " + vol[5] + "\n")
    f.write("ITEM: ATOMS id type x y z\n")
    
    for particle in pList:
        f.write(particle.get_LAMMPS() +  "\n")
        
    f.close()

def NewTCLScript(name,reactor):
    x,y,z = reactor.volume[1::2]
    
    f = open("sim_{}.tcl".format(name),"w+")
    f.write("set trajPath \"/home/stud/mb08/re83qin/projekt/Python/ParticleSimulation/\"\n")
    f.write("set simulationTraj \"LAMMPS_{}.lammpstrj\"\n".format(name))
    f.write("mol delete top \n")
    f.write("mol load lammpstrj $trajPath$simulationTraj \n")
    f.write("mol delrep 0 top \n")
    f.write("display resetview \n")
    
    for i in range(len(reactor.particles)):
        pType   = i + 1
        color   = i
        radius  = 0.7 * reactor.particles[i][2]
        
        f.write("mol representation VDW {} 16.000000 \n".format(radius)) #roughly 0.7*radius
        f.write("mol selection name {} \n".format(pType))
        #f.write("mol material Opaque \n")
        f.write("mol color ColorID {} \n".format(color))
        f.write("mol addrep top\n")
            
    f.write("set vol1 " + `x` + "\n")
    f.write("set vol2 " + `y` + "\n")
    f.write("set vol3 " + `z` + "\n")
    
    f.write("set minx -$vol1\n")
    f.write("set maxx $vol1\n")
    f.write("set miny -$vol2\n")
    f.write("set maxy $vol2\n")
    f.write("set minz -$vol3\n")
    f.write("set maxz $vol3\n")
    
    f.write("draw materials off\n")
    f.write("draw color white\n")
    
    f.write("draw line \"$minx $miny $minz\" \"$maxx $miny $minz\"\n")
    f.write("draw line \"$minx $miny $minz\" \"$minx $maxy $minz\"\n")
    f.write("draw line \"$minx $miny $minz\" \"$minx $miny $maxz\"\n")
    f.write("draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n")
    f.write("draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\"\n")
    f.write("draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\"\n")
    f.write("draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\"\n")
    f.write("draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\"\n")
    f.write("draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\"\n")
    f.write("draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\"\n")
    f.write("draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\"\n")
    f.write("draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\"\n")
    
    f.write("animate goto 0 \n")
    f.write("color Display Background white \n")
    f.close()


def NewMSDFile(name):
    open('MSD_{}.dat'.format(name),'w+')


def DumpOneMSDStep(name,reactor,time):
    pList = reactor.pList
    move = 0.0
    
    f = open('MSD_{}.dat'.format(name),'a')

    for particle in pList:
        init    = particle.initial
        loc     = particle.location
        for i in range(3):
            move += (init[i] - loc[i]) ** 2
    MSD = move / float(reactor.n)
    f.write(`time` +"  " + `MSD` +"  " +  `6*time*particle.diffusion` + "\n")
    

def NewMSDScript (name):
    f = open("msd.gnuplot","w+")
    f.write("set xlabel \"time[s]\"\n")
    f.write("set ylabel \"MSD\"\n")
    f.write("plot \"MSD_{}.dat\" using 1:2 with lines lw 5 , \\\n".format(name))
    f.write(" \"MSD_{}.dat\" using 1:3 with lines lw 5  \n".format(name))
    f.write("pause -1 \n")
    f.close()