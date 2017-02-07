import ROL
from firedrake import *

class FiredrakeLA(ROL.CustomLA):
    def __init__(self, vec, inner):
        ROL.CustomLA.__init__(self)
        self.data = vec
        self.inner = inner

    def plus(self, xx):
        self.data += xx.data

    def scale(self, alpha):
        self.data *= alpha

    def dot(self, xx):
        return self.inner.eval(self.data, xx.data)

    def clone(self):
        res = FiredrakeLA(self.data.copy(), self.inner)
        return res

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
y = Function(V)
z = Function(V)
uvec = u.vector()
yvec = y.vector()
zvec = z.vector()
x = FiredrakeLA(uvec, None)
y = FiredrakeLA(yvec, None)
z = FiredrakeLA(zvec, None)
x.checkVector(y, z)
