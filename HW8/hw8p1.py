import numpy as np
import matplotlib.pyplot as plt

# Utility functions

# Find the sum of values in an array arr taken to some power
def sumPow(arr, power):
    total = 0
    for i in arr:
        total += pow(i,power)
    return total

# Actual script

# Set up the random data to be fit
ns = 100
x = -1 + 5*np.random.random(ns)
y = 2*x*x*x - 4*x*x - x + 3 + 8*(-1 + 2*np.random.random(ns))

plt.plot(x,y,'o')
legendElements =  ["Data Points"]

maxDegree = 5;
for deg in range(1, maxDegree):
    
    # Define which degree to create a fit for
    degree = deg
    p = degree + 1
    a_rows, a_cols = (p,p)
    
    # See attached work for derivation
    # Need to define A matrix and b vector
    # This section creates the A matrix
    A = []
    for i in range(a_rows):
        col = []
        for j in range(a_cols):
            col.append(sumPow(x, p - j - 1 + i))
        A.append(col)
    
    # Now create the b vector
    b = []
    for i in range(p):
        total = 0
        for index in range(len(y)):
            total += y[index]*pow(x[index], i)
        b.append(total)
        
    # Now solve the system
    coeffs = np.linalg.solve(A, b)
    print('Coefficients for order ' , degree, ': ', coeffs)
    
    bestFitLine = np.zeros(len(y))
    for index in range(len(x)):
        value = 0;
        for coeff_index in range(len(coeffs)):
            value += coeffs[coeff_index]*pow(x[index], p - 1 - coeff_index)
        bestFitLine[index] = value
        
    # Since x values are unsorted, the plot looks bad. Need to sort
    # best fit y values by ascending x values.
    pairs = dict(zip(x, bestFitLine))
    sortedKeys = sorted(pairs)
    sortedBestFitLine = [pairs[keys] for keys in sortedKeys]
    plt.plot(sortedKeys, sortedBestFitLine)
    legendElements.append("Polynomial Degree " + str(degree))
    
    #nums = [1,5,2,7,3]
    #nums.sort()
    #print(nums)

plt.legend(legendElements)