import dolfin

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

import numpy as np
class dolfin_BasedLA(ROL.CustomLA):
    def __init__(self, n):
        ROL.CustomLA.__init__(self)
        self.data = dolfin.Vector(dolfin.mpi_comm_world(), n)
        self.clones = []

    def plus(self, x):
        self.data += x.data

    def scale(self, alpha):
        self.data *= alpha

    def __getitem__(self, i):
        return self.data[i][0]

    def __setitem__(self, i, v):
        self.data[i] = v

    def dot(self, x):
        val = self.data.inner(x.data)
        return val

    def clone(self):
        res = dolfin_BasedLA(self.data.size())
        res.data = self.data.copy()
        return res

x = dolfin_BasedLA(2)
y = dolfin_BasedLA(2)
z = dolfin_BasedLA(2)
x.data[0] = 1.0
x.data[1] = 1.5

x.checkVector(y, z)

algo.run(x, obj)
print x.data.array()
