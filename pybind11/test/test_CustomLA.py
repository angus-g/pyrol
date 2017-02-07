#!/usr/bin/env py.test
import ROL
import numpy as np


class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]

class MyObj2(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return -(x[0]+x[1] + x[0]*x[1])

    def gradient(self, g, x, tol):
        g[0] = -(1 + x[1]);
        g[1] = -(1 + x[0]);


class EqConstraint(ROL.EqualityConstraint):
    def __init__(self):
        ROL.EqualityConstraint.__init__(self)

    def value(self, c, x, tol):
        c[0] = x[0]*x[0]+x[1]*x[1]-1

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = 2*x[0]*v[0]+2*x[1]*v[1]

    def applyAdjointJacobian(self, ajv, v, x, tol):
        ajv[0] = v[0] * 2 * x[0]
        ajv[1] = v[0] * 2 * x[1]


class NPBasedLA(ROL.CustomLA):
    def __init__(self, size):
        ROL.CustomLA.__init__(self)
        self.data = np.zeros(size)
        self.size = size

    def plus(self, xx):
        self.data += xx.data

    def scale(self, alpha):
        self.data *= alpha

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, v):
        self.data[i] = v

    def dot(self, xx):
        return np.inner(self.data, xx.data)

    def dimension(self):
        return self.size

    def basis(self, i):
        res = NPBasedLA(self.size)
        res[i] = 1.0
        return res

    def clone(self):
        res = NPBasedLA(self.size)
        res.data = np.copy(self.data)
        return res

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

def test_checkVectorPassed():
    x = NPBasedLA(2)
    x.data[0] = 0.5
    x.data[1] = 0.5
    y = NPBasedLA(2)
    z = NPBasedLA(2)

    u = x.checkVector(y, z)
    assert sum(u) < 1e-12

def test_unconstrained():
    obj = MyObj()
    algo = ROL.Algorithm("Line Search", params)
    x = NPBasedLA(2)
    algo.run(x, obj)
    assert round(x[0] - 1.0, 8) == 0.0
    assert round(x[1], 8) == 0.0

def runBoundConstrained(algo):
    obj = MyObj()
    algo = ROL.Algorithm(algo, params)
    x = NPBasedLA(2)
    x_lo = NPBasedLA(2)
    x_lo[0] = -1
    x_lo[1] = -1
    x_up = NPBasedLA(2)
    x_up[0] = +0.7
    x_up[1] = +0.7
    bnd = ROL.BoundConstraint(x_lo, x_up, 1.0)
    algo.run(x, obj, bnd)
    assert round(x[0] - 0.7, 8) == 0.0
    assert round(x[1], 8) == 0.0

def test_boundConstrained():
    runBoundConstrained("Line Search")
    runBoundConstrained("Primal Dual Active Set")
    runBoundConstrained("Trust Region")

def test_equalityConstrainedCS():
    obj = MyObj2()
    algo = ROL.Algorithm("Composite Step", params)
    x = NPBasedLA(2)
    x[0] = 0.5 * 0.5**2
    x[1] = 0.5 * 0.5**2
    l = NPBasedLA(1)
    con = EqConstraint()
    algo.run(x, l, obj, con)
    assert round(x[0] - 0.707106, 5) == 0.0
    assert round(x[1] - 0.707106, 5) == 0.0

def test_equalityConstrainedAL():
    obj = MyObj2()
    algo = ROL.Algorithm("AugmentedLagrangian", params)
    x = NPBasedLA(2)
    x[0] = 0.5 * 0.5**2
    x[1] = 0.5 * 0.5**2
    con = EqConstraint()
    l = NPBasedLA(1)
    c = NPBasedLA(1)
    augLag = ROL.AugmentedLagrangian(obj, con, l, 10, x, c, params)
    bnd = ROL.BoundConstraint() # inactive bound, here because ROL needs it
    algo.run(x, l, augLag, con, bnd)
    assert round(x[0] - 0.707106, 5) == 0.0
    assert round(x[1] - 0.707106, 5) == 0.0
