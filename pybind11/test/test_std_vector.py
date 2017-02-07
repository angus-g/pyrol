#!/usr/bin/env py.test

import ROL

# Create a quadratic objective function around (1.0, 0.0)
class Objective(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]
obj = Objective()

def test_std_vector():

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
    x = ROL.StdVector(2)
    x[0] = -1.0
    x[1] = 2.0

    algo.run(x, obj)
    # Check answer to 8 decimal places
    assert round(x[0] - 1.0, 8) == 0.0
    assert round(x[1], 8) == 0.0
