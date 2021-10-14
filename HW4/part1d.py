import numpy as np
from time import process_time #because my python is 3.8

N = 100000000
r = np.arange(N, dtype = np.double) # make sure data is stored as a double
t1 = process_time()

sum = 0.0
for i in range(N):
    sum += r[i]*r[i]

mag = np.sqrt(sum)
t2 = process_time()
print("mag = %g"%mag)
print("Loop took %g seconds"%(t2 - t1))

#Numpy vector version

t1 = process_time()
mag = np.sqrt(np.dot(r,r))
t2 = process_time()
print("NumPy vectors took %g seconds"%(t2 - t1))
