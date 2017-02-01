#!/usr/bin/python
import ROL

print dir(ROL)

class GreatObjective(ROL.ObjectiveDouble):
    def __init__(self):
        pass
    def value(self, vec, tol):
        print vec, tol
        return vec.norm()

obj = GreatObjective()
v = ROL.StdVectorDouble([1, 2, 3, 5])
print obj, dir(obj)
print obj.value(v, 1e-6)

vdot = ROL.StdVectorDouble([0, 0, 0, 0])

# print obj.gradient(v, vdot, 0.01)


print v.norm()
print v.dot(v)

algo = ROL.AlgorithmDouble("hello", {'param':True})

algo.run(v, obj)
