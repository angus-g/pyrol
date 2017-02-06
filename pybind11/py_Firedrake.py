import ROL
print "import good"
from firedrake import *

class FiredrakeLA(ROL.CustomLA):
    def __init__(self, vec, inner):
        self.data = vec
        self.inner = inner

    def plus(self, xx):
        # pass
        self.data += xx.data

    def scale(self, alpha):
        # print "here"
        self.data *= alpha

    def dot(self, xx):
        return self.inner.eval(self.data, xx.data)

    def clone(self):
        res = FiredrakeLA(self.data.copy(), self.inner)
        return res

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
uvec = u.vector()
x = FiredrakeLA(uvec, None)
