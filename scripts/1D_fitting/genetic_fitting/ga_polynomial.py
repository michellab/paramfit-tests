
from pyevolve import G1DList, GSimpleGA, Selectors, Scaling, DBAdapters
from random import seed, randint, random

def eval_polynomial(x, *coefficients):
    result = 0
    for exponent, coeff in enumerate(coefficients):
     result += coeff*x**exponent

    return result

def generate_fitness_function(sample_points):
    def fitness_function(chromosome):
         score = 0
         for point in sample_points:
             print(point[0])
             delta = abs(eval_polynomial(point[0], *chromosome) - point[1])
             score += delta
         score = -score
         return score
    return fitness_function


# Generate a random polynomial, and generate sample points from it
seed()

source_polynomial = []
for i in xrange(randint(1, 5)):
 source_polynomial.append(randint(-20,20))

sample_points = []
for i in xrange(20):
 n = randint(-100, 100)
 #*souce_polynomial are the coeficient of the polynoms

 sample_points.append((n, eval_polynomial(n, *source_polynomial)))
 #n = x points   poly = y points
# Create the population
genome = G1DList.G1DList(5)
genome.evaluator.set(generate_fitness_function(sample_points))
genome.setParams(rangemin=-50, rangemax=50)

# Set up the engine
ga = GSimpleGA.GSimpleGA(genome)
ga.setPopulationSize(1000)
ga.selector.set(Selectors.GRouletteWheel)

# Change the scaling method
pop = ga.getPopulation()
pop.scaleMethod.set(Scaling.SigmaTruncScaling)

# Start the algorithm, and print the results.
ga.evolve(freq_stats=10)
#print(ga.bestIndividual())
best = ga.bestIndividual()
print("Source polynomial: " + repr(source_polynomial))
print("Genetic points:")
print(best[0])
print(best[1])
print(best[2])
print(best[3])
print(best[4])

#print("Sample points: " + repr(sample_points))
