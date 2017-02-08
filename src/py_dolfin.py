#!/usr/bin/env python
from math import cos, pi
import dolfin
import ROL

class Objective(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)
        self.n = 6.0

    def value(self, x, tol):
        val = 0.0
        n = self.n
        for i in range(x.data.local_size()):
            val += (x[i] - i)**n
        return val

    def gradient(self, g, x, tol):
        n = self.n
        for i in range(x.data.local_size()):
            g[i] = n*(x[i] - i)**(n-1)

obj = Objective()

params = """
<ParameterList>
  <ParameterList name="Step">
    <ParameterList name="Line Search">
      <ParameterList name="Descent Method">
        <Parameter name="Type" type="string" value="Quasi-Newton Method"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  <ParameterList name="Status Test">
    <Parameter name="Gradient Tolerance" type="double" value="1e-20"/>
    <Parameter name="Step Tolerance" type="double" value="1e-20"/>
    <Parameter name="Iteration Limit" type="int" value="1000"/>
  </ParameterList>
</ParameterList>
"""

algo = ROL.Algorithm("Line Search", params)

class dolfin_BasedLA(ROL.CustomLA):
    def __init__(self, n):
        ROL.CustomLA.__init__(self)
        self.data = dolfin.Vector(dolfin.mpi_comm_world(), n)
        self.size = n

    def plus(self, x):
        self.data += x.data

    def scale(self, alpha):
        self.data *= alpha

    def __getitem__(self, i):
        return self.data[i][0]

    def __setitem__(self, i, v):
        self.data[i] = v

    def dot(self, x):
        return self.data.inner(x.data)

    def dimension(self):
        return self.data.size()

    def dimension(self):
        return self.size

    def basis(self, i):
        res = dolfin_BasedLA(self.size)
        res[i] = 1.0
        return res

    def clone(self):
        return dolfin_BasedLA(self.data.size())
#        res = dolfin_BasedLA(self.data.size())
#        res.data = self.data.copy()
#        return res

x = dolfin_BasedLA(20)
y = dolfin_BasedLA(20)
z = dolfin_BasedLA(20)
x.data[0] = 1.0
x.data[1] = 1.5

x.checkVector(y, z)

algo.run(x, obj)
import numpy
numpy.set_printoptions(precision=4, suppress=True)
print x.data.array()
