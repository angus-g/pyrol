#!/usr/bin/python
import ROL
print dir(ROL)

v = ROL.StdVector(20)
print v.dimension()

class GreatObjective(ROL.Objective):
    def __init__(self):
        pass
    def value(self, x, tol):
        print x, tol
        return x.norm()

w = GreatObjective()
print w.value(v, 0.1)

algo = ROL.Algorithm("work", {'example':'param'})
print algo
