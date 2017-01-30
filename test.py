#!/usr/bin/python
import ROL
print dir(ROL.MyVector)

v = ROL.MyVector(20)

obj = ROL.MyObjective()
print v

print obj.valuezz(v, 0.1555)

print v.dimension()
