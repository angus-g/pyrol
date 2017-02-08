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
x = ROL.StdVector(2)
x[0] = -1.0
x[1] = 2.0

print x[0]
print x[1]
algo.run(x, obj)
print x[0]
print x[1]
