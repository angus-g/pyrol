from firedrake import UnitSquareMesh, FunctionSpace, Function, \
    Expression, TrialFunction, TestFunction, inner, grad, DirichletBC, \
    dx, Constant, solve, assemble, as_backend_type, File
from firedrake_LA import FiredrakeLA
import ROL

n = 64
mesh = UnitSquareMesh(n, n)
use_correct_riesz_and_inner = True
outdir = "output_riesz_%s/" % use_correct_riesz_and_inner
V = FunctionSpace(mesh, "Lagrange", 1)  # space for state variable
M = FunctionSpace(mesh, "DG", 0)  # space for control variable
beta = 1e-4
yd = Function(V)
# yd.interpolate(Expression("(x[0] <= 0.5)*(x[1] <= 0.5)", degree=1))
yd.interpolate(Expression("(1.0/(2*pi*pi)) * sin(pi*x[0]) * sin(pi*x[1])",
               degree=3))


def solve_state(u):
    y = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(y), grad(v)) * dx
    L = u * v * dx
    bc = DirichletBC(V, Constant(0.0), "on_boundary")
    y = Function(V)
    solve(a == L, y, bc)
    return y


def solve_adjoint(u, y):
    lam = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(lam), grad(v)) * dx
    L = -(y-yd) * v * dx
    bc = DirichletBC(V, Constant(0.0), "on_boundary")
    lam = Function(V)
    solve(a == L, lam, bc)
    return lam


class L2Inner(object):

    def __init__(self):
        self.A = assemble(TrialFunction(M)*TestFunction(M)*dx)
        self.Ap = as_backend_type(self.A).mat()
        self.bcs = [DirichletBC(M, Constant(0.0), "on_boundary")]

    def eval(self, _u, _v):
        upet = as_backend_type(_u).vec()
        vpet = as_backend_type(_v).vec()
        A_u = self.Ap.createVecLeft()
        self.Ap.mult(upet, A_u)
        return vpet.dot(A_u)

    def riesz_map(self, derivative):
        rhs = Function(M, val=derivative.dat)
        res = Function(M)
        solve(self.A, res, rhs, bcs=self.bcs,
              solver_parameters={
                  'ksp_monitor': False,
                  'ksp_rtol': 1e-9, 'ksp_atol': 1e-10, 'ksp_stol': 1e-16,
                  'ksp_type': 'cg', 'pc_type': 'hypre',
                  'pc_hypre_type': 'boomeramg'
              })
        return res.vector()


state_file = File(outdir + "state.pvd")
control_file = File(outdir + "control.pvd")


class Objective(ROL.Objective):
    '''Subclass of ROL.Objective to define value and gradient for problem'''
    def __init__(self, inner_product):
        ROL.Objective.__init__(self)
        self.inner_product = inner_product
        self.u = Function(M)
        self.y = Function(V)

    def value(self, x, tol):
        u = self.u
        y = self.y
        return assemble(0.5 * (y-yd) * (y-yd) * dx + 0.5 * beta * u * u * dx)

    def gradient(self, g, x, tol):
        u = self.u
        y = self.y
        lam = solve_adjoint(u, y)
        v = TestFunction(M)
        L = beta * u * v * dx - lam * v * dx
        deriv = assemble(L)
        grad = self.inner_product.riesz_map(deriv)
        grad.dat.copy(g.vec.dat)

    def update(self, x, flag, iteration):
        u = Function(M, x.vec)
        u.vector().dat.copy(self.u.vector().dat)
        control_file.write(self.u)
        y = solve_state(self.u)
        y.vector().dat.copy(self.y.vector().dat)
        state_file.write(self.y)


parametersXML = """
<ParameterList>
  <ParameterList name="Step">
    <ParameterList name="Line Search">
      <ParameterList name="Descent Method">
        <Parameter name="Type" type="string"
          value="Quasi-Newton Method"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  <ParameterList name="Status Test">
    <Parameter name="Gradient Tolerance" type="double" value="1e-8"/>
    <Parameter name="Step Tolerance" type="double" value="1e-8"/>
    <Parameter name="Iteration Limit" type="int" value="4"/>
  </ParameterList>
</ParameterList>
"""
params = ROL.ParameterList(parametersXML)
inner_product = L2Inner()
obj = Objective(inner_product)
u = Function(M)
opt = FiredrakeLA(u.vector(), inner_product)
algo = ROL.Algorithm("Line Search", params)
algo.run(opt, obj)
File("res.pvd").write(Function(M, opt.vec))
