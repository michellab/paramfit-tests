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


def sample_points_generation(angles):

    sampled = []
    for angle in angles:
        y = 1. + 3.*math.cos(angle)
        x = angle
        sampled.append((x,y))

    return sampled

def eval_function(x, *coefficients):

    result = 0
    #print(coefficients)
    for a_val, b_val in enumerate(coefficients):
     result += a_val + b_val*math.cos(x)

    return result

def generate_fitness_function(sample_points):
    def fitness_function(chromosome):
         score = 0
         for point in sample_points:
             delta = (eval_function(point[0], *chromosome) - point[1])**2
             score += delta
         score = -score
         return score
    return fitness_function


def comparison(a_val,b_val):
    #This will give me the comparison beween fitted and real function
    #betahat list (2 elems) is the final results we  obtained from GA
    y_values = []
    y_fit    = []
    x_data_points = np.linspace(0,2*math.pi,100)
    for x in x_data_points:
        y = 1+3*math.cos(x)
        y_values.append(y)

        fitting = a_val + b_val*math.cos(x)
        y_fit.append(fitting)

    plt.figure(1)
    xx = x_data_points
    plt.scatter(xx,y_values,color="b",label="Original")
    plt.scatter(xx,y_fit,color="r",label="fitting")
    plt.legend(loc="best")
    plt.show()

#########
#MAIN###
#####
random.seed()
n_input = 20
angles = np.linspace(0,2*math.pi,n_input)
sample_points = sample_points_generation(angles)


#initialize alleles
setOfAlleles=GAllele.GAlleles()

#first variable
a_val=np.arange(0,10,0.1)
setOfAlleles.add(GAllele.GAlleleList(a_val))
#second variable
b_val=np.arange(0,10,0.1)
setOfAlleles.add(GAllele.GAlleleList(b_val))
print(setOfAlleles)
#initialize genome with defined alleles
genome = G1DList.G1DList(2)
genome.setParams(allele=setOfAlleles)

#define evaluator function
genome.evaluator.set(generate_fitness_function(sample_points))
genome.mutator.set(Mutators.G1DListMutatorAllele)
genome.initializator.set(Initializators.G1DListInitializatorAllele)

ga = GSimpleGA.GSimpleGA(genome)
ga.setPopulationSize(500)
ga.selector.set(Selectors.GRouletteWheel)
ga.setGenerations(1000)
ga.setCrossoverRate(0.20)
ga.setElitism(True)
# Start the algorithm, and print the results.
# Change the scaling method
pop = ga.getPopulation()
pop.scaleMethod.set(Scaling.ExponentialScaling)
ga.evolve(freq_stats=1000)

best=ga.bestIndividual()


print(best)
