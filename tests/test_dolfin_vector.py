#!/usr/bin/env python
import dolfin
import ROL
from ROL.dolfin_vector import DolfinVector

class Objective(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)
        self.n = 4.0

    def value(self, x, tol):
        val = 0.0
        n = self.n
        for i in range(x.vec.local_size()):
            val += (x[i] - i)**n
        return val

    def gradient(self, g, x, tol):
        n = self.n
        for i in range(x.vec.local_size()):
            g[i] = n*(x[i] - i)**(n-1)

obj = Objective()


paramsdict = {
    "Step": {
        "Line Search": {
            "Descent Method": {
                "Type": "Quasi-Newton Method"
            }
        }
    },
    "Status Test":{
        "Gradient Tolerance" :"1e-20",
        "Step Tolerance" : "1e-20",
        "Iteration Limit": "100",
    }
}
params = ROL.ParameterList(paramsdict, "Parameters")

algo = ROL.Algorithm("Line Search", params)

def test_dolfin_la_checkvector():
    x = DolfinVector(dolfin.Vector(dolfin.MPI.comm_world, 2))
    y = DolfinVector(dolfin.Vector(dolfin.MPI.comm_world, 2))
    z = DolfinVector(dolfin.Vector(dolfin.MPI.comm_world, 2))
    x.vec[0] = 1.0
    x.vec[1] = 1.5

    x.checkVector(y, z)

def test_dolfin_la_objective():
    n = 10
    m = n*(n-1)/2
    x = DolfinVector(dolfin.Vector(dolfin.MPI.comm_world, n))
    algo.run(x, obj)
    assert abs((sum(x.vec.get_local()) - m)/m) < 1e-6
