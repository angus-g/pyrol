from dolfin import *
from ROLUtils.dolfin_LA import dolfinLA as LA
import numpy
import ROL

E = Constant(1.0) # Young's modulus
nu = Constant(0.3) # Poisson's ratio
lmbda = (E*nu) / ((1+nu)*(1-2*nu)) # formulae for lmbda, mu in terms of E, nu taken from http://en.wikipedia.org/wiki/Lam%C3%A9_parameters
mu = E / (2*(1 + nu))

V = 0.3 # Maximum volume allowed
p = Constant(4.0)
gamma = Constant(1e-6)
beta = Constant(1e-8)

def alpha(rho):
    return rho**p + gamma

def alphadash(rho):
    return p * (rho**(p-1))

def eps(v): return sym(grad(v))

N = 100
mesh = RectangleMesh(Point(0.0, 0.0), Point(1.6, 1.0), N, N)
V0 = FunctionSpace(mesh, "CG", 1) # for the control
V1 = VectorFunctionSpace(mesh, "CG", 1) # for the displacement

class SourceExpression(Expression):
    def eval(self, values, x):
        values[0] = 0.0
        values[1] = 0.0

        if all(x == (1.6, 0.5)):
            values[1] = -1.0

    def value_shape(self):
        return (2,)

f = interpolate(SourceExpression(degree=1), V1)

class InitialExpressionA(Expression):
    def eval(self, values, x):
        values[0] = V

def solve_state(rho):
    """Solve the forward problem for a given fluid distribution rho(x)."""
    u = TrialFunction(V1)
    v = TestFunction(V1)

    F = alpha(rho) * 2 * mu * inner(eps(u), eps(v)) * dx
    F += lmbda * div(u) * div(v) * dx
    F -= inner(f, v) * dx
    bcs = [DirichletBC(V1, (0.0, 0.0), "x[0] < DOLFIN_EPS")]
    u = Function(V1)
    solve(lhs(F) == rhs(F), u, bcs=bcs)
    return u

def solve_adjoint(rho, u):
    v = TrialFunction(V1)
    du = TestFunction(V1)
    F = alpha(rho) * 2 * mu * inner(eps(du), eps(v)) * dx
    F += lmbda * div(du) * div(v) * dx
    F += inner(f, du) * dx
    bcs = [DirichletBC(V1, (0.0, 0.0), "x[0] < DOLFIN_EPS")]
    adj = Function(V1)
    solve(lhs(F) == rhs(F), adj, bcs=bcs)
    return adj

class L2Inner(object):

    def __init__(self):
        self.A = assemble(TrialFunction(V0)*TestFunction(V0)*dx)

    def eval(self, _u, _v):
        A_u = _v.copy()
        self.A.mult(_u, A_u)
        return _v.inner(A_u)

    def riesz_map(self, derivative):
        res = Function(V0)
        rhs = Function(V0, derivative)
        solve(self.A, res.vector(), rhs.vector())
        return res.vector()

dot_product = L2Inner()

state_file = File("state.pvd")
control_file = File("control.pvd")
class ObjR(ROL.Objective):
    '''Subclass of ROL.Objective to define value and gradient for problem'''
    def __init__(self, inner_product):
        ROL.Objective.__init__(self)
        self.inner_product = inner_product
        self.rho = Function(V0)
        self.state = Function(V1)

    def value(self, x, tol):
        rho = self.rho
        u = self.state
        return assemble(inner(f, u) * dx + beta * inner(grad(u), grad(u)) * dx)

    def gradient(self, g, x, tol):
        rho = self.rho
        u = self.state
        v = solve_adjoint(rho, u)
        drho = TestFunction(V0)
        L = alphadash(rho) * drho * 2 * mu * inner(eps(u), eps(v)) * dx + 2 * beta * inner(grad(rho), grad(drho)) * dx
        deriv = assemble(L)
        if self.inner_product is not None:
            gra = self.inner_product.riesz_map(deriv)
        else:
            gra = deriv
        g.scale(0)
        g.vec += gra

    def update(self, x, flag, iteration):
        rho = Function(V0, x.vec)
        self.rho.assign(rho)
        state = solve_state(self.rho)
        self.state.assign(state)
        control_file << self.rho
        state_file << self.state

class VolConstraint(ROL.EqualityConstraint):

    def __init__(self, inner_product):
        ROL.EqualityConstraint.__init__(self)
        self.inner_product = inner_product

    def value(self, cvec, xvec, tol):
        rho = Function(V0, xvec.vec)
        val = assemble(rho * dx) - V
        cvec[0] = val

    def applyJacobian(self, jv, v, x, tol):
        drho = Function(V0, v.vec)
        jv[0] = assemble(drho * dx)

    def applyAdjointJacobian(self, ajv, v, x, tol):
        drho = TestFunction(V0)
        deriv = assemble(drho*dx)
        if self.inner_product is not None:
            gra = self.inner_product.riesz_map(deriv)
        else:
            gra = deriv
        ajv.scale(0)
        ajv.vec += gra
        ajv.scale(v[0])

# Initialise 'ROLVector'
l = ROL.StdVector(1)
c = ROL.StdVector(1)
v = ROL.StdVector(1)
v[0] = 1.0
dualv = ROL.StdVector(1)
v.checkVector(c, l)

x = interpolate(InitialExpressionA(degree=1), V0)
x = LA(x.vector(), dot_product)
g = Function(V0)
g = LA(g.vector(), dot_product)
d = interpolate(Expression("1 + x[0] * (1-x[0])*x[1] * (1-x[1])", degree=1), V0)
d = LA(d.vector(), dot_product)
x.checkVector(d, g)

jd = Function(V0)
jd = LA(jd.vector(), dot_product)

lower = interpolate(Constant(0.0), V0)
lower = LA(lower.vector(), dot_product)
upper = interpolate(Constant(1.0), V0)
upper = LA(upper.vector(), dot_product)

# Instantiate Objective class for poisson problem
obj = ObjR(dot_product)
obj.checkGradient(x, d, 7, 1)
volConstr = VolConstraint(dot_product)
volConstr.checkApplyJacobian(x, d, jd, 6, 1)
volConstr.checkAdjointConsistencyJacobian(v, d, x)

with open('input.xml', 'r') as myfile:
    parametersXML=myfile.read().replace('\n', '')
set_log_level(30)

params = ROL.ParameterList(parametersXML)
bound_constraint = ROL.BoundConstraint(lower, upper, 1.0)

alg2 = ROL.Algorithm("Augmented Lagrangian", params)
penaltyParam = 1
augLag = ROL.AugmentedLagrangian(obj, volConstr, l, penaltyParam, x, c, params)
alg2.run(x, l, augLag, volConstr, bound_constraint)
