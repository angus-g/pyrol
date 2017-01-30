#!/usr/bin/python
import ROL

v = ROL.MyVector(200)
print v, dir(v)

print v.dimension()

obj = ROL.MyObjective()
print obj, dir(obj)

print obj.value(v, 0.1555)
