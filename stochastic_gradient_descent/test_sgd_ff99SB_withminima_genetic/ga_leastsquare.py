#Jan 2017 Stefano Bosisio
#Here I want to fit the  function  y(x) = 1. + 3.*cos(x) by using a GA
#to minimize LLS

import numpy as np
import math, random
import matplotlib.pyplot as plt
from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Mutators
from pyevolve import Initializators
from pyevolve import Scaling
from pyevolve import GAllele



def eval_function(eqm,emm,phi_angle,psi_angle,phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi):
# eval : energies sum_i^m(EQM - EMM - ct - all the dihedrals to fit)^2

    sum_squared = (eqm - emm  - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos - k3psi*math.cos(3*psi_angle) - k3psi*psi_3cos)**2

    #print(sum_cos)
    return sum_squared

def generate_fitness_function(dict_structures,phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos):
    def fitness_function(chromosome):
        #print("Function")
        score = 0
        for key in dict_structures:
            eqm = dict_structures[key][0]
            emm = dict_structures[key][1]
            phi_angle = dict_structures[key][2]
            psi_angle = dict_structures[key][3]
            delta = eval_function(eqm,emm,phi_angle,psi_angle,phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos, *chromosome)
            score += delta
        print(score)
        score = -score
        return score
    return fitness_function



#########
#MAIN###
#####
random.seed()

#read the file and recreate th edictionary, we will use it for genetic purposes
input_f = open("all_data.dat","r").readlines()
#dict_structures[n_struc] = [eqm,emm,phi,psi]
#eqm in kcal/mol and relative, so the same for emm, phi and psi already in rad
#emm already with the offset
dict_structures = {}
stop = 43
for lines in input_f :
    #now here hard cod eplease or I get tired
    splitter = lines.split()
    n_struct = int(splitter[0])
    eqm = float(splitter[2])
    emm = float(splitter[3])
    phi = float(splitter[4])
    psi = float(splitter[5])
    dict_structures[n_struct] = [eqm,emm,phi,psi]
    if n_struct==stop:
        break
#print(dict_structures)
phi_1cos = 0.570087
phi_2cos = 0.205221
phi_3cos = -0.845488
psi_1cos = 0.221390
psi_2cos = -0.011197
psi_3cos = -0.433653

k1_phi = -0.033287
k2_phi = 0.367068
k3_phi = 0.329716

k1_psi = -0.29444
k2_psi = -1.569341
k3_psi = -0.583119


#initialize alleles
setOfAlleles=GAllele.GAlleles()

#k1_phi
end = k1_phi + (k1_phi*0.5)  #a 10% modification
start = k1_phi - (k1_phi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen1_phi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen1_phi))
#k2_phi
end = k2_phi + (k2_phi*0.5)  #a 10% modification
start = k2_phi - (k2_phi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen2_phi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen2_phi))

#k3_phi
end = k3_phi + (k3_phi*0.5)  #a 10% modification
start = k3_phi - (k3_phi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen3_phi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen3_phi))

#k1_psi
end = k1_psi + (k1_psi*0.5)  #a 10% modification
start = k1_psi - (k1_psi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen1_psi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen1_psi))

#k2_psi
end = k2_psi + (k2_psi*0.5)  #a 10% modification
start = k2_psi - (k2_psi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen2_psi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen2_psi))

#k3_psi
end = k3_psi + (k3_psi*0.5)  #a 10% modification
start = k3_psi - (k3_psi*0.5)
if start > end:
    tmp = end
    end = start
    start =  tmp
gen3_psi=np.arange(start,end,0.0001)
setOfAlleles.add(GAllele.GAlleleList(gen3_psi))


#initialize genome with defined alleles
genome = G1DList.G1DList(6)
genome.setParams(allele=setOfAlleles)

#define evaluator function
genome.evaluator.set(generate_fitness_function(dict_structures,phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos))
genome.mutator.set(Mutators.G1DListMutatorAllele)
genome.initializator.set(Initializators.G1DListInitializatorAllele)

ga = GSimpleGA.GSimpleGA(genome)
ga.setPopulationSize(500)
ga.selector.set(Selectors.GRouletteWheel)
ga.setGenerations(1000)
ga.setCrossoverRate(0.2)
ga.setElitism(True)
# Start the algorithm, and print the results.
# Change the scaling method
pop = ga.getPopulation()
pop.scaleMethod.set(Scaling.ExponentialScaling)
ga.evolve(freq_stats=1000)

best=ga.bestIndividual()


print(best)
