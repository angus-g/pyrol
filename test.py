#!/usr/bin/python
import ROL
print dir(ROL.MyVector)

v = ROL.MyVector(20)

obj = ROL.MyObjective()
print obj.value(v)

print v.dimension()
