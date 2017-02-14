from dolfin import *
from dolfin_LA import dolfinLA as LA
import numpy
import ROL

mu = Constant(1.0)                   # viscosity
alphaunderbar = 2.5 * mu / (100**2)  # parameter for \alpha
alphabar = 2.5 * mu / (0.01**2)      # parameter for \alpha
q = Constant(0.5) # q value that controls difficulty/discrete-valuedness of solution

def alpha(rho):
    """Inverse permeability as a function of rho, equation (40)"""
    return alphabar + (alphaunderbar - alphabar) * rho * (1 + q) / (rho + q)

def alphadash(rho):
    return (alphaunderbar - alphabar) * (1 * (1 + q) / (rho + q) - rho * (1 + q)/((rho + q)*(rho + q)))

N = 100
delta = 1.5  # The aspect ratio of the domain, 1 high and \delta wide
V = (1.0/3) * delta  # want the fluid to occupy 1/3 of the domain

mesh = RectangleMesh(mpi_comm_world(), Point(0.0, 0.0), Point(delta, 1.0), N, N)
A = FunctionSpace(mesh, "DG", 0)        # control function space

U_h = VectorElement("CG", mesh.ufl_cell(), 2)
P_h = FiniteElement("CG", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, U_h*P_h)          # mixed Taylor-Hood function space

class InflowOutflow(Expression):
    def eval(self, values, x):
        values[1] = 0.0
        values[0] = 0.0
        l = 1.0/6.0
        gbar = 1.0

        if x[0] == 0.0 or x[0] == delta:
            if (1.0/4 - l/2) < x[1] < (1.0/4 + l/2):
                t = x[1] - 1.0/4
                values[0] = gbar*(1 - (2*t/l)**2)
            if (3.0/4 - l/2) < x[1] < (3.0/4 + l/2):
                t = x[1] - 3.0/4
                values[0] = gbar*(1 - (2*t/l)**2)

    def value_shape(self):
        return (2,)

def solve_state(rho):
    """Solve the forward problem for a given fluid distribution rho(x)."""
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)

    F = (alpha(rho) * inner(u, v) * dx + inner(grad(u), grad(v)) * dx +
         inner(grad(p), v) * dx  + inner(div(u), q) * dx + Constant(0) * q * dx)
    bc = DirichletBC(W.sub(0), InflowOutflow(degree=1), "on_boundary")
    w = Function(W)
    solve(lhs(F) == rhs(F), w, bcs=bc)
    return w

def solve_adjoint(rho, w): 
    (u, p) = split(w)
    (v, q) = TrialFunctions(W)
    (du, dp) = TestFunctions(W)
    F = inner(alpha(rho) * u, du) * dx
    F += inner(2 * mu * grad(u), grad(du)) * dx 
    F += inner(alpha(rho) * du, v) * dx 
    F += inner(mu * grad(du), grad(v)) * dx 
    F += inner(grad(dp), v) * dx 
    F += inner(div(du), q) * dx
    bc = DirichletBC(W.sub(0), Constant((0,0)), "on_boundary")
    adj = Function(W)
    solve(lhs(F) == rhs(F), adj, bcs=bc)
    return adj

class L2Inner(object):

    def __init__(self):
        self.A = assemble(TrialFunction(A)*TestFunction(A)*dx)
        self.Ap = as_backend_type(self.A).mat()

    def eval(self, _u, _v):
        upet = as_backend_type(_u).vec()
        vpet = as_backend_type(_v).vec()
        A_u = self.Ap.createVecLeft()
        self.Ap.mult(upet, A_u)
        return vpet.dot(A_u)

    def riesz_map(self, derivative):
        res = Function(A)
        rhs = Function(A, derivative)
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
        self.rho = Function(A)
        self.state = Function(W)

    def value(self, x, tol):
        rho = self.rho
        state = self.state
        (u, p) = split(state)
        return assemble(0.5 * inner(alpha(rho) * u, u) * dx + mu * inner(grad(u), grad(u)) * dx)

    def gradient(self, g, x, tol):
        rho = self.rho
        state = self.state
        (u, p) = split(state)
        lam = solve_adjoint(rho, state)
        (v, q)= split(lam)
        drho = TestFunction(A)
        L = 0.5 * alphadash(rho) * drho * inner(u, u) * dx + alphadash(rho) * drho * inner(u, v) * dx
        deriv = assemble(L)
        if self.inner_product is not None:
            grad = self.inner_product.riesz_map(deriv)
        else:
            grad = deriv
        g.scale(0)
        g.vec += grad

    def update(self, x, flag, iteration):
        rho = Function(A, x.vec)
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
        a = Function(A, xvec.vec)
        val = assemble(a * dx) - V
        cvec[0] = val

    def applyJacobian(self, jv, v, x, tol):
        da = Function(A, v.vec)
        jv[0] = assemble(da * dx)

    def applyAdjointJacobian(self, ajv, v, x, tol):
        da = TestFunction(A)
        deriv = assemble(da*dx)
        if self.inner_product is not None:
            grad = self.inner_product.riesz_map(deriv)
        else:
            grad = deriv
        ajv.scale(0)
        ajv.vec += grad
        ajv.scale(v[0])

# Initialise 'ROLVector'
l = ROL.StdVector(1)
c = ROL.StdVector(1)
v = ROL.StdVector(1)
v[0] = 1.0
dualv = ROL.StdVector(1)
v.checkVector(c, l)

x = interpolate(Constant(V/delta), A)
x = LA(x.vector(), dot_product)
g = Function(A)
g = LA(g.vector(), dot_product)
d = interpolate(Expression("1 + x[0] * (1-x[0])*x[1] * (1-x[1])", degree=1), A)
d = LA(d.vector(), dot_product)
x.checkVector(d, g)

jd = Function(A)
jd = LA(jd.vector(), dot_product)

lower = interpolate(Constant(0.0), A)
lower = LA(lower.vector(), dot_product)
upper = interpolate(Constant(1.0), A)
upper = LA(upper.vector(), dot_product)

# Instantiate Objective class for poisson problem
obj = ObjR(dot_product)
# obj.checkGradient(x, d, 4, 2)
volConstr = VolConstraint(dot_product)
volConstr.checkApplyJacobian(x, d, jd, 3, 1)
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
