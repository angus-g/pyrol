#!/usr/bin/env python
import dolfin
import ROL 
class Objective(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)
        self.n = 4.0

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

parameterXML = """
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
    <Parameter name="Iteration Limit" type="int" value="100"/>
  </ParameterList>
</ParameterList>
"""
params = ROL.ParameterList(parameterXML)

algo = ROL.Algorithm("Line Search", params)

class dolfinLA(ROL.CustomLA):
    def __init__(self, vec):
        ROL.CustomLA.__init__(self)
        self.data = vec 

    def plus(self, x):
        # PETSc doesn't like adding vector to themselves,
        # so we have to add this check
        xpet = dolfin.as_backend_type(x.vec).vec()
        thispet = dolfin.as_backend_type(self.vec).vec()
        if xpet == thispet:
            self.vec *= 2
        else:
            self.vec += x.vec

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
        res = dolfinLA(self.size)
        res[i] = 1.0
        return res

    def clone(self):
        return dolfinLA(self.data.copy())

def test_dolfin_la_checkvector():
    x = dolfinLA(dolfin.Vector(dolfin.mpi_comm_world(), 2))
    y = dolfinLA(dolfin.Vector(dolfin.mpi_comm_world(), 2))
    z = dolfinLA(dolfin.Vector(dolfin.mpi_comm_world(), 2))
    x.data[0] = 1.0
    x.data[1] = 1.5
    
    x.checkVector(y, z)

def test_dolfin_la_objective():
    n = 10
    m = n*(n-1)/2
    x = dolfinLA(dolfin.Vector(dolfin.mpi_comm_world(), n))
    algo.run(x, obj)
    assert abs((sum(x.data.array()) - m)/m) < 1e-6
