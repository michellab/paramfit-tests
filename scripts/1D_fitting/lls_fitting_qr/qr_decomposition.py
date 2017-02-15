import numpy as np
import math, random
import matplotlib.pyplot as plt


def random_points(number):

    x_data_points = []
    y_data_points = []
    y_data_squared= []

    for i in range(0,number):
        ran_x = random.uniform(0,2*math.pi)
        func  = 1. + (3.*math.cos(ran_x))  # the function I want to fit is cosine
        x_data_points.append(ran_x)
        y_data_points.append(func)
        product = func*func
        y_data_squared.append(product)

    #plt.figure(1)
    #xx = x_data_points
    #plt.scatter(xx,y_data_points,color="b")
    #plt.show()


    return x_data_points,y_data_points, y_data_squared

def comparison(betahat):

    y_values = []
    y_fit    = []
    x_data_points = np.linspace(0,2*math.pi,100)
    for x in x_data_points:
        y = 1+3*math.cos(x)
        y_values.append(y)

        fitting = betahat[0] + betahat[1]*math.cos(x)
        y_fit.append(fitting)

    plt.figure(1)
    xx = x_data_points
    plt.scatter(xx,y_values,color="b",label="Original")
    plt.scatter(xx,y_fit,color="r",label="fitting")
    plt.legend(loc="best")
    plt.show()



n_input = 20
x_data,y_data, y_squared = random_points(n_input)

elem_a = n_input + 1
elem_b = sum(y_data)
elem_c = sum(y_data)
elem_d = 0
for val in x_data:
    elem_d += math.cos(val)**2

matrix_list = [elem_a,elem_b,elem_c,elem_d]
matrix = np.array(matrix_list).reshape(2, 2)
qr = np.linalg.qr(matrix,mode="complete")
print("Q matrix\n")
q = qr[0]
print(qr[0])
print("R matrix\n")
r = qr[1]
print(qr[1])

sum_y = sum(y_data)
sum_y_cos = 0
for i,val in enumerate(x_data,0):
    sum_y_cos += y_data[i]*math.cos(val)

y_list = [sum_y,sum_y_cos]
y = np.array(y_list,).reshape(-1, 1)
betaHat = np.linalg.inv(r).dot(q.T).dot(y)

print("b values\n")
print(betaHat)
#now test on S :
S = 0.0
for i,val in enumerate(x_data,0):
    S += (y_data[i] - (betaHat[0] + betaHat[1]*math.cos(val)))**2

if
#comparison(betaHat)
