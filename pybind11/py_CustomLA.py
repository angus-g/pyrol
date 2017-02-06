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
class NPBasedLA(ROL.CustomLA):
    def __init__(self):
        ROL.CustomLA.__init__(self)
        self.data = np.asarray([0.0,0.0])
        self.clones = []

    def plus(self, xx):
        self.data[0] += xx.data[0]
        self.data[1] += xx.data[1]

    def scale(self, alpha):
        self.data[0] *= alpha
        self.data[1] *= alpha

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, v):
        self.data[i] = v

    def dot(self, xx):
        return self.data[0]*xx.data[0]+self.data[1]*xx.data[1]

    def clone(self):
        res = NPBasedLA()
        res.data[0] = self.data[0]
        res.data[1] = self.data[1]
        return res

x = NPBasedLA()
y = NPBasedLA()
z = NPBasedLA()
x.data[0] = 1.0
x.data[1] = 1.5

x.checkVector(y, z)

algo.run(x, obj)
print x.data
