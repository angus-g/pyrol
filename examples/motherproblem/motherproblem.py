# from firedrake import UnitSquareMesh, FunctionSpace, Function, \
    # Expression, TrialFunction, TestFunction, inner, grad, DirichletBC, \
    # dx, Constant, solve, assemble, as_backend_type, File
# from firedrake_LA import FiredrakeLA as LA
# backend = "firedrake"
from dolfin_LA import dolfinLA as LA
from dolfin import UnitSquareMesh, FunctionSpace, Function, \
    Expression, TrialFunction, TestFunction, inner, grad, DirichletBC, \
    dx, Constant, solve, assemble, as_backend_type, File, CellFunction, refine
backend = "dolfin"

import ROL

n = 16
mesh = UnitSquareMesh(n, n)

if backend == "dolfin":
    import numpy
    numpy.random.seed(0)
    def randomly_refine(initial_mesh, ratio_to_refine= .3):
       cf = CellFunction('bool', initial_mesh)
       for k in range(len(cf)):
           if numpy.random.rand() < ratio_to_refine:
               cf[k] = True
       return refine(initial_mesh, cell_markers = cf)

    k = 2
    for i in range(k):
       mesh = randomly_refine(mesh)
else:
    k = 0

use_correct_riesz_and_inner = True

outdir = "output_riesz_%s_refinements_%i/" % (use_correct_riesz_and_inner, k)

V = FunctionSpace(mesh, "Lagrange", 1)  # space for state variable
M = FunctionSpace(mesh, "DG", 0)  # space for control variable
beta = 1e-4
yd = Function(V)
yd.interpolate(Expression("(1.0/(2*pi*pi)) * sin(pi*x[0]) * sin(pi*x[1])",
               degree=3))
# uncomment this for a 'more difficult' target distribution
# yd.interpolate(Expression("(x[0] <= 0.5)*(x[1] <= 0.5)", degree=1))


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
        self.bcs = [DirichletBC(M, Constant(0.0), "on_boundary")]

    def eval(self, _u, _v):
        y = _v.copy()
        y.zero()
        self.A.mult(_u, y)
        return _v.inner(y)

    def riesz_map(self, derivative):
        if backend == "firedrake":
            rhs = Function(M, val=derivative.dat)
            res = Function(M)
            solve(self.A, res, rhs, bcs=self.bcs)
            # solve(self.A, res, rhs, bcs=self.bcs,
            #       solver_parameters={
            #           'ksp_monitor': False,
            #           'ksp_rtol': 1e-9, 'ksp_atol': 1e-10, 'ksp_stol': 1e-16,
            #           'ksp_type': 'cg', 'pc_type': 'hypre',
            #           'pc_hypre_type': 'boomeramg'
            #       })
            return res.vector()
        else:
            self.bcs[0].apply(self.A)
            res = Function(M)
            rhs = Function(M, derivative)
            solve(self.A, res.vector(), rhs.vector())

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
        if self.inner_product is not None:
            grad = self.inner_product.riesz_map(deriv)
        else:
            grad = deriv
        g.scale(0)
        g.vec += grad

    def update(self, x, flag, iteration):
        u = Function(M, x.vec)
        self.u.assign(u)
        y = solve_state(self.u)
        self.y.assign(y)
        if backend == "firedrake":
            control_file.write(self.u)
            state_file.write(self.y)
        else:
            control_file << self.u
            state_file << self.y


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
if use_correct_riesz_and_inner:
    inner_product = L2Inner()
else:
    inner_product = None
obj = Objective(inner_product)
u = Function(M)
opt = LA(u.vector(), inner_product)
xlo = Function(M)
xlo.interpolate(Constant(0.0))
x_lo = LA(xlo.vector(), inner_product)
xup = Function(M)
xup.interpolate(Constant(0.9))
x_up = LA(xup.vector(), inner_product)
bnd = ROL.BoundConstraint(x_lo, x_up, 1.0)

algo = ROL.Algorithm("Line Search", params)
algo.run(opt, obj, bnd)
if backend == "firedrake":
    File("res.pvd").write(Function(M, opt.vec))
else:
    File("res.pvd") << Function(M, opt.vec)
