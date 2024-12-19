import numpy as np

a = np.arange(12).reshape((4,3))

print(a)

ix = np.ix_([0,0],[2,2])

print(a[ix])


