from virSim import *

maxBirthProb = 0.9
clearProb = 0.05
popDensity = 0.1
maxPop = 1000
mutProb = 0.005

def test_getResistPop():
    num_viruses = 10
    resistances = {'guttagonol':True}
    viruses = []
    drugResist = ['guttagonol']
    for i in range(num_viruses):
        virus = ResistantVirus(maxBirthProb,clearProb,resistances,mutProb)
        viruses.append(virus)

    patientX = Patient(viruses,maxPop)
    resVir = patientX.getResistPop(drugResist)
    return resVir
    
#print test_getResistPop()

def test_reproduce():
    resistances = {'guttagonol':True,'grimpex':False}
#    drugs = ['guttagonol','grimpex']
#    drugs = []
    drugs = ['grimpex'] # if virus is not resistant to any drug in drugs then it does not reproduce
    num_viruses = 5
    viruses = []
    for i in range(num_viruses):
        virus = ResistantVirus(maxBirthProb,clearProb,resistances,mutProb)
        print 'created virus ', virus
        viruses.append(virus)
    vir_pop = viruses[:]
    for virus in viruses:
        print '\n'
        print 'virus, ',virus
        try:
            
            offspring = virus.reproduce(popDensity,drugs )
            vir_pop.append(offspring)
          
        except NoChildException:
            continue
    return len(vir_pop)
    
#print test_reproduce()



def test_update():
    num_viruses = 10
    print 'num_viruses ', num_viruses, 
    viruses = []
    population = []
    resistances = {'guttagonol':True}
    numSteps = 5
    print 'numSteps ',numSteps
    
    for i in range(num_viruses):
        virus = ResistantVirus(maxBirthProb,clearProb,resistances,mutProb)
        viruses.append(virus)
        print 'virus ',i , 'name ', virus
        print 'self.resistance ', virus.resistances
    patientX = Patient(viruses,maxPop)
    patientX.addPrescription('guttagonol')
    
    
    for i in range(1,numSteps+1):
         pop = patientX.update()
         population.append(pop)

    return population

#print test_update()

