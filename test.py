#!/usr/bin/python
import ROL

print dir(ROL)

v = ROL.MyVector(200)
print v, dir(v)

print v.dimension()

obj = ROL.MyObjective()
print obj, dir(obj)

print obj.value(v, 0.1555)

v = ROL.StdVectorDouble([1,2])
print v.norm()
print v.dot(v)
