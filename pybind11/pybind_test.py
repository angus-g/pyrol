import ROL
class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2
    
    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]
obj = MyObj()
algo = ROL.Algorithm("Line Search", {"test":"test"})

x = ROL.StdVector(2)
x[0] = -1.0
x[1] = 2.0

print x[0] 
print x[1] 
algo.run(x, obj)
print x[0] 
print x[1] 


_x = ROL.EigenVector(2)
_x[0] = -1.0
_x[1] = 2.0

print _x[0] 
print _x[1] 
algo = ROL.Algorithm("Line Search", {"test":"test"})
algo.run(_x, obj)
print _x[0] 
print _x[1] 

import numpy as np
a = np.asarray([0, 1, 2, 3])
evec = ROL.EigenVector(a)
print evec[0]
print evec[1]
print evec[2]
print evec[3]
a[3] = 5.0
print evec[3]
evec[3]=3.1
print a[3]

