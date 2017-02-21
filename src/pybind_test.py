#!/usr/bin/env python

import ROL
class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]
obj = MyObj()

params = """
<ParameterList>
  <ParameterList name="Step">
    <ParameterList name="Line Search">
      <ParameterList name="Descent Method">
        <Parameter name="Type" type="string" value="Quasi-Newton Method"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  <ParameterList name="Status Test">
    <Parameter name="Gradient Tolerance" type="double" value="1e-12"/>
    <Parameter name="Step Tolerance" type="double" value="1e-16"/>
    <Parameter name="Iteration Limit" type="int" value="10"/>
  </ParameterList>
</ParameterList>
"""

algo = ROL.Algorithm("Line Search", params)

import numpy as np
def iplus(x, y):
    x.data[0] += y.data[0]
    x.data[1] += y.data[1]

def iscale(x, alpha):
    x.data[0] *= alpha

class NPBasedLA(ROL.CustomLA):
    def __init__(self):
        ROL.CustomLA.__init__(self)
        # super(ROL.CustomLA, self).__init__(self)
        self.data = np.asarray([0, 0])
        self.clones = []
        print "construct"

    def iplus(self, xx):
        print "here plus"

    def iscale(self, alpha):
        print "here scale"

    def iclone(self):
        print "here clone"
        res = NPBasedLA()
        self.clones.append(res)
        # res.iscale(1.0)
        return res

x = NPBasedLA()
y = NPBasedLA()
z = NPBasedLA()
# x.scale(3.0)
# x.plus(y)
# # xx = x.clone()
# # xx.scale(1.0)
# yy = x.iclone()
# yy.iscale(1.0)
x.checkClone()
x.checkVector(y, z)
# # algo.run(x, obj)
