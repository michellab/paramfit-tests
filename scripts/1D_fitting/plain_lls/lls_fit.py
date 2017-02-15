import numpy as np
import matplotlib.pyplot as plt
import math,random

def reference_function(x_points) :
    fx = []
    for x in x_points:
        val = math.cos(x)
        fx.append(val)
    return fx

def random_datapoints(fx):

    n_data = 10
    x_data_points = []
    y_data_points = []

    for  i in range(0,n_data):
        valx = random.uniform(0,2*math.pi)
        x_data_points.append(valx) # x
        valy = math.cos(valx)
        y_data_points.append(valy)

    return x_data_points, y_data_points

def plot_fitting(fx,betaHat,x_points,input_datapoints):

    plt.figure(1)
    xx = x_points
    yy = np.array(betaHat[0] + betaHat[1] * xx)
    plt.plot(xx, yy.T, color='b')
    plt.scatter(input_datapoints[:, 0], input_datapoints[:, 1], color='r')
    plt.show()


#main#
x_axis = np.linspace(0,2*math.pi,100)
fx = reference_function(x_axis)
x_data,y_data = random_datapoints(fx)
array_values = []
for i, val in enumerate(x_data,0):
    array_values.append([val,y_data[i]])

xy_array = np.array(array_values)
m = len(xy_array)

X = np.array([np.ones(m), xy_array[:, 0]]).T

y = np.array(xy_array[:, 1]).reshape(-1, 1)


betaHat = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)
print("Beta hat solution\n")
print(betaHat)

plot_fitting(fx,betaHat,x_axis,xy_array)
