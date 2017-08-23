import numpy
import random
import pylab

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. 
    """    


class SimpleVirus(object):
    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        
        clearProb: Maximum clearance probability (a float between 0-1).
        """
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb
        
        
        
    def doesClear(self):
        """
        Stochastically determines whether this virus is cleared from the
        patient's body at a time step. 

        returns: Using a random number generator (random.random()), this method
        returns True with probability self.clearProb and otherwise returns
        False.
        """
        random_float = random.random()
        if random_float <= self.clearProb:
            return True
        else:
            return False
        
    
    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the SimplePatient and
        Patient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """

        probability = self.maxBirthProb*(1 - popDensity)
        random_float = random.random()
        if  random_float <= probability:
            return SimpleVirus(self.maxBirthProb, self.clearProb)
        else:
            raise NoChildException(Exception)
            
class SimplePatient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """
    
    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)
        
        maxPop: the  maximum virus population for this patient (an integer)
        """

        self.viruses = viruses
        self.maxPop = maxPop
        

    def getTotalPop(self):
        """
        Gets the current total virus population. 

        returns: The total virus population (an integer)
        """
     
        return len(self.viruses)

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:

        - Determine whether each virus particle survives and updates the list
          of virus particles accordingly.

        - The current population density is calculated. This population density
          value is used until the next call to update() 

        - Determine whether each virus particle should reproduce and add
          offspring virus particles to the list of viruses in this patient.                    

        returns: the total virus population at the end of the update (an
        integer)
        """
        virus_pop = self.viruses[:]
        for virus in self.viruses:
            if virus.doesClear():
                virus_pop.remove(virus)
        popDensity = len(virus_pop)/float(self.maxPop )      
        self.viruses = virus_pop[:] 
        for virus in virus_pop:
            try:
                offspring = virus.reproduce(popDensity)
                self.viruses.append(offspring)
            except NoChildException:
                continue
        return self.getTotalPop()    
        
def problem2():
    """
    Run the simulation and plot the graph for problem 2 (no drugs are used,
    viruses do not have any drug resistance).    

    Instantiates a patient, runs a simulation for 300 timesteps, and plots the
    total virus population as a function of time.    
    """

    maxBirthProb = 0.1
    clearProb = 0.05
    viruses = []
    maxPop = 1000
    time = 200
    timesteps = pylab.arange(1,time+1)
    for i in range(0,100):
        virus=SimpleVirus(maxBirthProb,clearProb)
        viruses.append(virus)

       
    patient = SimplePatient(viruses,maxPop)
    population = []
    for i in timesteps:
        current_pop = patient.update()
        population.append(current_pop)
    pylab.plot(timesteps,population)
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population')
    pylab.title('Virus population as a function of time')
    pylab.show()

#print problem2()
    

class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """    
    
    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.
        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        
        clearProb: Maximum clearance probability (a float between 0-1).
        
        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'grimpex',False}, means that this virus
        particle is resistant to neither guttagonol nor grimpex.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.        
        """
  
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb
        self.resistances = resistances
        self.mutProb = mutProb
        
    def getResistance(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in Patient to determine how many virus
        particles have resistance to a drug.        

        drug: the drug (a string).

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
  
        if self.resistances[drug]:
            return True
        else:
            return False

          
    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient class.

        If the virus particle is not resistant to any drug in activeDrugs,
        then it does not reproduce. Otherwise, the virus particle reproduces
        with probability:       
        
        self.maxBirthProb * (1 - popDensity).                       
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). 

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.        

        For example, if a virus particle is resistant to guttagonol but not
        grimpex, and `self.mutProb` is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90% 
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        grimpex and a 90% chance that the offspring will not be resistant to
        grimpex.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population        

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings). 
        
        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.         
        """
#        print 'reproduce function'
        resistance = True

        for drug in activeDrugs:
 #           print 'DRUG: ',drug
            if not self.getResistance(drug):
#                print 'drug ',drug
                resistance = False
#                print 'resistance ',resistance
##            else:
##                resistance = False

#        print 'resistance ',resistance

        if resistance or len(activeDrugs)==0:
            probability = self.maxBirthProb*(1-popDensity)
            randomFloat = random.random()
            if randomFloat <= probability:
                randomRes = random.random()
#                print 'new offspring'
                if randomRes<=self.mutProb:
##                    print 'resistance CHANGED'
                    offres = {}
                    for resistance in self.resistances:
                        if self.resistances[resistance]:
                            offres[resistance]=False
                        else:
                            offres[resistance]=True
                    return ResistantVirus(self.maxBirthProb,self.clearProb,offres,self.mutProb)
##                else:
##                    print 'resistance unchanged'
                return ResistantVirus(self.maxBirthProb, self.clearProb,self.resistances,self.mutProb)
            else:
#                print 'Exception error: Not Born'
                raise NoChildException(Exception)            
        else:
#            print 'Exception error: No resistance to Drug'
            raise NoChildException(Exception)
        
        
            
class Patient(SimplePatient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """
    
    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).               

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)
        
        maxPop: the  maximum virus population for this patient (an integer)
        """

        self.viruses = viruses
        self.maxPop = maxPop
        self.drugs = []
        
    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the 
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: list of drugs being administered to a patient is updated
        """

        if newDrug not in self.drugs:
            self.drugs.append(newDrug)
        print 'self.drugs ',self.drugs
        return self.drugs

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """

        return self.drugs

        
    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in 
        drugResist.        

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'grimpex'])

        returns: the population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """

        resistPop = []
        for virus in self.viruses:
            resistance = False
            for drug in drugResist:
                if virus.getResistance(drug):
                    resistance = True
                else:
                    resistance = False
            if resistance:
                resistPop.append(virus)
                
        return len(resistPop)
    
    
    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of 
          virus particles accordingly
          
        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Determine whether each virus particle should reproduce and add
          offspring virus particles to the list of viruses in this patient. 
          The listof drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces. 

        returns: the total virus population at the end of the update (an
        integer)
        """
 
##        print '\n'
##        print 'update function'
        virus_pop = self.viruses[:]

        for virus in self.viruses:
            if virus.doesClear():
                virus_pop.remove(virus)
        popDensity = len(virus_pop)/float(self.maxPop )      
        self.viruses = virus_pop[:]
        count = 1
        for virus in virus_pop:
##            print '\n'
##            print 'virus',virus
            try:
##                print '\n'
##                print 'try ',count
                count+=1
                offspring = virus.reproduce(popDensity,self.drugs)
                self.viruses.append(offspring)
            except NoChildException:
                continue
            
        return self.getTotalPop()    
                

def problem4():
    """
    Runs simulations and plots graphs for problem 4.

    Instantiates a patient, runs a simulation for 150 timesteps, adds
    guttagonol, and runs the simulation for an additional 150 timesteps.

    total virus population vs. time  and guttagonol-resistant virus population
    vs. time are plotted
    """

    maxPop = 1000
    virNum = 100
    viruses = []
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = {'guttagonol':False}
    mutProb = 0.005
    print 'Creating a virus database'
    for i in range(virNum):
        virus = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)
        viruses.append(virus)

    print 'Initializating a Patient.'
    patientX = Patient(viruses,maxPop)
    time = 150
    timeSteps = pylab.arange(1,time+1)
    population = []
    print 'Running simulation for',time,'timesteps.'
    for j in timeSteps:
        current_pop = patientX.update()
#        print current_pop
        population.append(current_pop)

    print 'Creating a plot for ',time,' timesteps.'
    pylab.plot(population)
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population')
    pylab.title('Total virus population vs. time')
    pylab.show()    

    
    print 'Administering guttagonol to the patient.'
    patientX.addPrescription('guttagonol')

    print 'Running simulation for',time,'timesteps.'
    for j in timeSteps:
        current_pop = patientX.update()
#        print current_pop
        population.append(current_pop)

    print 'Creating a plot for ',time,' timesteps.'  
    pylab.plot(population)
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population')
    pylab.title('Guttagonol-resistant virus population vs. time')   
    pylab.show()

##    print 'Administering grimpex to the patient.'
##    patientX.addPrescription('grimpex')
##
##    print 'Running simulation for',time,'timesteps.'
##    for j in timeSteps:
##        current_pop = patientX.update()
###        print current_pop
##        population.append(current_pop)
##
##    print 'Creating a plot for ',time,' timesteps.'  
##    pylab.plot(population)
##    pylab.xlabel('Time Steps')
##    pylab.ylabel('Virus Population')
##    pylab.title('Grimpex-resistant virus population vs. time')   
##    pylab.show()    
        
#problem4()

      
def problem5():
    """
    Runs simulations and make histograms for problem 5.

    Runs multiple simulations to show the relationship between delayed treatment
    and patient outcome.

    Histograms of final total virus populations are displayed for delays of 300,
    150, 75, 0 timesteps (followed by an additional 150 timesteps of
    simulation).    
    """
    steps = [300,150,75,0]
    addSteps = 150
    maxPop = 1000
    virNum = 10
    viruses = []
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = {'guttagonol':False}
    mutProb = 0.005
    numPatients = 10
    print 'Creating a virus database'
    for i in range(virNum):
        virus = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)
        viruses.append(virus)

#    finalPop = []
    for i in steps:
        print 'Running simulation for ',i,' timesteps.'
        result = []
        for j in range(numPatients):
#            print 'sim ',j
            patientX = Patient(viruses,maxPop)
            timeSteps = pylab.arange(1,i+1)
            population = []
            for k in timeSteps:
                current_pop = patientX.update()
                population.append(current_pop)
            patientX.addPrescription('guttagonol')
            for h in range(addSteps):
                current_pop = patientX.update()
                population.append(current_pop)
            totalPop = patientX.getTotalPop()
            result.append(totalPop)

        print 'Printing Virus population for',numPatients,'patients.'
        print result
#        avgPop = sum(result)/float(len(result))
#        print avgPop
#        finalPop.append(avgPop)
        pylab.figure()
        print 'Creating a histogram for ',i,' timesteps.'
        title = 'A final virus population with a drug after first '+ str(i)+' time steps'
        pylab.hist(result)
        pylab.xlabel('Total Virus Population')
        pylab.ylabel('Number of patients')
        pylab.title(title)
        pylab.show()
        

#problem5()        

def problem6():
    """
    Runs simulations and make histograms for problem 6.

    Runs multiple simulations to show the relationship between administration
    of multiple drugs and patient outcome.
    
    Histograms of final total virus populations are displayed for lag times of
    150, 75, 0 timesteps between adding drugs (followed by an additional 150
    timesteps of simulation).
    """

    virNum = 100
    maxPop = 1000
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = {'guttagonol':False,'grimpex':False}
    mutProb = 0.005
    viruses = []
    virPop = []
    simSteps = [300,150,75,0]
    addSteps = 150
    numPatients = 30
    
    print 'Creating a virus database.'
    for i in range(virNum):
        virus = ResistantVirus(maxBirthProb,clearProb, resistances, mutProb)
        viruses.append(virus)
        
    print 'Initializating a Patient.'
    patientX = Patient(viruses,maxPop)
    
    firstTimeSteps = 150
    print 'Running simulation for first ',firstTimeSteps,' timesteps.'
    for i in range(firstTimeSteps):
        current_pop = patientX.update()
        virPop.append(current_pop)
    print virPop

    print 'Administering guttagonol to the patient.'
    patientX.addPrescription('guttagonol')
#    virPop.append('Gutt')

    
    for j in simSteps:
        print 'Running simulation for following condition : ',j,'timesteps.'
        result = virPop[:]
        finalPop = []

        for k in range(numPatients):
#            print 'Creating Patient Y'
            population = result[:]
            patientY = Patient(patientX.viruses,maxPop)
#            print '--running simulation for ',j,' timesteps.'
            for m in range(j):
                pop = patientY.update()
                population.append(pop)
            
            patientY.addPrescription('grimpex')
#            population.append('Grim')
            for n in range(addSteps):
                pop = patientY.update()
                population.append(pop)
#            print 'final population'
#            print population
            totalPop = patientY.getTotalPop()
            finalPop.append(totalPop)
        
        print 'Printing Virus population for',numPatients,'patients.'
        print finalPop
        pylab.figure()
        print 'Creating a histogram for ',j,' timesteps.'
        title = 'A final virus population with second drug after '+ str(j)+' time steps ater fist drug.'
        pylab.hist(finalPop)
        pylab.xlabel('Total Virus Population')
        pylab.ylabel('Number of patients')
        pylab.title(title)
        pylab.show()
      
#problem6()


     
def problem7(guttagonolAt,grimpexAt):
    """
    Run simulations and plot graphs examining the relationship between
    administration of multiple drugs and patient outcome.

    Plots of total and drug-resistant viruses vs. time are made for a
    simulation with a 300 time step delay between administering the 2 drugs and
    a simulations for which drugs are administered simultaneously.        
    """

    virNum = 100
    maxPop = 1000
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = {'guttagonol':False,'grimpex':False}
    mutProb = 0.005
    viruses = []
    totViruses = []
    guttagonolResistant = []
    grimpexResistant = []
    dualResistant = []
    timeStep = 1
#    guttagonolAt = 150
#    grimpexAt = 450 
    
    print 'Creating a virus database.'
    for i in range(virNum):
        virus = ResistantVirus(maxBirthProb,clearProb, resistances, mutProb)
        viruses.append(virus)
        
    print 'Initializating a Patient.'
    patientX = Patient(viruses,maxPop)

    print 'Running simulation'
    while timeStep <= grimpexAt + 150:
        virPop = patientX.update()
        timeStep += 1        
        if timeStep == guttagonolAt:
            patientX.addPrescription('guttagonol')
        if timeStep == grimpexAt:
            patientX.addPrescription('grimpex')
        totViruses.append(virPop)
        guttagonolResistant.append(patientX.getResistPop(['guttagonol']))
        grimpexResistant.append(patientX.getResistPop(['grimpex']))
        dualResistant.append(patientX.getResistPop(['guttagonol','grimpex']))

    print 'Creating a plot of Total Population.'

    pylab.figure()
    pylab.plot(totViruses, label='Total')
    pylab.plot(guttagonolResistant, label= 'Guttagonol resistant virus')
    pylab.plot(grimpexResistant, label= 'Grimpex resistant virus')
    pylab.plot(dualResistant, label= 'Resistant to both drugs')
    pylab.xlabel('Time Steps')
    pylab.ylabel('Virus Population')
    pylab.title('Total virus population vs. time')
    pylab.legend(loc = 4)
    pylab.show()

#problem7(guttagonolAt = 150, grimpexAt = 450)
problem7(guttagonolAt = 150, grimpexAt = 150)

