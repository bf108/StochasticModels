#  Simulating the Spread of Disease and Virus Population Dynamics 

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

    def getMaxBirthProb(self):
        """
        Returns the max birth probability.
        """
        return self.maxBirthProb

    def getClearProb(self):
        """
        Returns the clear probability.
        """

        return self.clearProb

    def doesClear(self):
        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """

        # Produce random float between 0 - 1.
        # Compare float to getClearProb. If < TRUE, Else False
        if self.getClearProb() >= random.random():
            return True
        else:
            return False

    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
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

        if self.maxBirthProb * (1-popDensity) >= random.random():
            return SimpleVirus(self.getMaxBirthProb(),self.getClearProb())
        else:
            raise NoChildException(Exception)


class Patient(object):
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

        maxPop: the maximum virus population for this patient (an integer)
        """
   
        self.viruses = viruses
        self.maxPop = maxPop

    def getViruses(self):
        """
        Returns the viruses in this Patient.
        """
    
        return self.viruses

    def getMaxPop(self):
        """
        Returns the max population.
        """
        
        return self.maxPop

    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
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
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        
        surviving_virus = []
        for virus in self.getViruses():
            if virus.doesClear() == False:
                surviving_virus.append(virus)
        
        popDensity = len(surviving_virus)/float(self.getMaxPop())
        
        new_virus = []
        for virus in surviving_virus:
            try:
                new_virus.append(virus.reproduce(popDensity))
            except NoChildException:
                pass

                
        self.viruses = new_virus + surviving_virus
        return len(self.viruses)

def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph for problem 3 (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    
    Typical Values Used for Calculation:
    maxBirthProb: 0.1
    clearProb: 0.05
        
    """

    steps = 301
    virus = []
    for _ in range(numViruses):
        virus.append(SimpleVirus(maxBirthProb,clearProb))

    patient = Patient(virus,maxPop)
    
    print("Created Viruses")
    
    time = range(0,steps,1)
    results = {}
    for i in range(numTrials):
        virus_count = []
        # print('Trial ', i)
        for _ in range(steps):
            patient.update()
            virus_count.append(patient.getTotalPop())
            # len(virus_count)
        results[i] = virus_count
    mean = []
    for i in range(steps):
        val = 0
        for j in range(numTrials):
            val += results[j][i]     
        mean.append(val/float(numTrials))
    pylab.plot(time,mean,'r',label = "No of Virus")
    pylab.title('SimpleVirus simulation')
    pylab.xlabel("Time Steps")
    pylab.ylabel("Average Virus Population")
    pylab.legend(loc = 'best')
    pylab.show()


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
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """
        
        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        self.resistance = resistances
        self.mutProb = mutProb

    def getResistances(self):
        """
        Returns the resistances for this virus.
        """
        
        return self.resistance

    def getMutProb(self):
        """
        Returns the mutation probability for this virus.
        """
        
        return self.mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        return self.getResistances().get(drug)
        
    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.maxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """

        # Checks to see whether virus is resistant to all drugs listed
        resistant = True
        for drug in activeDrugs:
            if self.getResistances().get(drug) == False:
                resistant = False
                break
        # Only reproduces if virus is resistant to all drugs
        # and probabilty maxBirth * (1-popDensity)
        if resistant == True and self.maxBirthProb * (1-popDensity) >= random.random():
            #If the virus reproduces. This assigns the new resistance to offspring
            new_res = {}
            for drug in self.resistance.items():
                if drug[1] == True:
                    if (1-self.mutProb) > random.random():
                        new_res[drug[0]] = True
                    else:
                        new_res[drug[0]] = False
                else:
                    if self.mutProb > random.random():
                        new_res[drug[0]] = True
                    else:
                        new_res[drug[0]] = False
            return ResistantVirus(self.getMaxBirthProb(),self.getClearProb(),new_res,self.getMutProb())
        else:
            raise NoChildException(Exception)
            

class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """

        Patient.__init__(self,viruses,maxPop)
        self.presDrugs = []

    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """
        
        if newDrug not in self.presDrugs:
            self.presDrugs.append(newDrug)   

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """
        
        return self.presDrugs

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        
        count = 0
        for virus in self.getViruses():
            sub_count = 0
            for drug in drugResist:
                if virus.getResistances().get(drug) == True:
                    sub_count += 1
            if sub_count == len(drugResist):
                count += 1
    
        return count

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """
        
        surviving_virus = []
        for virus in self.getViruses():
            if virus.doesClear() == False:
                surviving_virus.append(virus)
        
        popDensity = len(surviving_virus)/float(self.getMaxPop())
        
        new_virus = []
        for virus in surviving_virus:
            try:
                new_virus.append(virus.reproduce(popDensity, self.getPrescriptions()))
            except NoChildException:
                pass
                
        self.viruses = new_virus + surviving_virus
        return len(self.viruses)


def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, resistances, mutProb, numTrials):
    """
    Runs simulations and plots graphs for problem 5.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})

    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
    """
    
    steps = 300
    virus = []
    for _ in range(numViruses):
        virus.append(ResistantVirus(maxBirthProb,clearProb,resistances,mutProb))
 
    print("Created Viruses")
    
    time = range(0,steps,1)
    results = {}
    results2 = {}
    
    for i in range(numTrials):
        patient = TreatedPatient(virus,maxPop)
        virus_count = []
        res_virus = []
        
        for j in range(steps):
            if j == 50:
                for item in resistances.keys():
                    patient.addPrescription(item)
            patient.update()
            virus_count.append(patient.getTotalPop())
            res_virus.append(patient.getResistPop([item for item in resistances.keys()]))
        results[i] = virus_count
        results2[i] = res_virus
    
    mean = []
    mean2 = []
    for i in range(steps):
        val = 0
        val2 = 0
        for j in range(numTrials):
            val += results[j][i]
            val2 += results2[j][i]
        mean.append(val/float(numTrials))
        mean2.append(val2/float(numTrials))

    pylab.plot(time,mean,'r',label = "No of Virus")
    pylab.plot(time,mean2,'b',label = "No of Virus Resistant to G")
    pylab.title('ResistantVirus simulation')
    pylab.xlabel("Time Steps")
    pylab.ylabel("# viruses")
    pylab.legend(loc = 'best')
    pylab.show()

''' 

    Test Simulation of Simple Virus with no drugs to combat virus population 
    
    This is done with the following parameters:
        Num Virsus to Start (Postive Int): 100
        Max Virus Population (Postive Int): 5000
        Max Birth Rate (Probability Virus will Reproduce) (float between 0-1 ): 0.05
        Clear Probabilty (Probability Virsu will Die) (float between 0-1): 0.03
        Number of Trials used to estimate mean outcome (Postive Int): 10
        
'''

simulationWithoutDrug(100,5000,0.05,0.03,10)


''' 

    Test Simulation of Drug Resitnt Virus and Antidotes to combat virus population 
    
    This is done with the following parameters:
        Num Virsus to Start (Postive Int): 100
        Max Virus Population (Postive Int): 5000
        Max Birth Rate (Probability Virus will Reproduce) (float between 0-1 ): 0.1
        Clear Probabilty (Probability Virsu will Die) (float between 0-1): 0.05
        Resistances (Dictionary of drugs which virus population 
        will be treated with and current resistance status) e.g 'penecilin':False
        Number of Trials used to estimate mean outcome (Postive Int): 10
        
'''

simulationWithDrug(100,5000,0.05,0.03,{'penecilin':False,'MMR':False},0.03,10)




