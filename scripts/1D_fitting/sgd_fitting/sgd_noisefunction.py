#JAN 2017 Stefano Bosisio
#Stochastic gradient descent fitting

import numpy as np
import math, random
import matplotlib.pyplot as plt

def generate_sample(number):

    x_axis = np.linspace(0.1,20,number)
    sampled = []

    for x in x_axis:
        funct = 1.473 + 3.842*math.sin(x) + 4.729*math.sin(4*x) + 5.980*math.log(5*x)
        sampled.append((x,funct))

    return sampled

def evaluate_function(a_val,b_val,c_val,d_val,sampled_points):

    #here evaluate sampled_points - newfunction **2
    sum_squared = 0
    for val in sampled_points:
        y_val = val[1]
        x_val = val[0]
        sum_squared += (y_val - a_val - b_val*math.sin(x_val) - c_val*math.sin(4*x_val) - d_val*math.log(5*x_val))**2

    return sum_squared
#MAIN#

n_input = 10 #number of input points
samples = generate_sample(n_input) #list of x and y points for the fitting
convergence = True #convergence will be reached if S<10^-9 ?
#print(samples)
A_t = 0 #starting point for A
B_t = 0 #starting point for B
C_t = 0
D_t = 0
alpha =0.01# 1./n_input # training parameter usually 0.1
print("Alpha %.4f\n" %alpha)
T_period = 100000
counter = 0
s_tmp = 1.
while (convergence):

    #select a random sampled point
    random_points = random.randint(0,n_input-1)
    x_val = samples[random_points][0]
    y_val = samples[random_points][1]
    #A(t) = A(t-1) - alpha*grad(lossfunction)
    #grad w.r.t A
    A_new = A_t - alpha*                   (A_t + B_t*math.sin(x_val) + C_t*math.sin(4*x_val) + D_t*math.log(5*x_val)- y_val)
    B_new = B_t - alpha*(math.sin(x_val)*  (A_t + B_t*math.sin(x_val) + C_t*math.sin(4*x_val) + D_t*math.log(5*x_val) - y_val))
    C_new = C_t - alpha*(math.sin(4*x_val)*(A_t + B_t*math.sin(x_val) + C_t*math.sin(4*x_val) + D_t*math.log(5*x_val) - y_val))
    D_new = D_t - alpha*(math.log(5*x_val)*(A_t + B_t*math.sin(x_val) + C_t*math.sin(4*x_val) + D_t*math.log(5*x_val) - y_val))

    s_val = evaluate_function(A_new,B_new,C_new,D_new,samples)

    if s_val < 0.0000001:
        print("This are the values found for fitting")
        print(A_new)
        print(B_new)
        print(C_new)
        print(D_new)
        convergence = False
    else:

        print(A_new,B_new,C_new,D_new)
        print(s_val)
        print(counter)
        A_t = A_new
        B_t = B_new
        C_t = C_new
        D_t = D_new
        counter+=1
