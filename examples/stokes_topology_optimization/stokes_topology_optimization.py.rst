Stokes Topology Optimization
============================

We consider the topology optimization example from `dolfin-adjoint <http://www.dolfin-adjoint.org/en/latest/documentation/stokes-topology/stokes-topology.html/>`_ ::

    from dolfin import *
    from ROL.dolfin_vector import DolfinVector as FeVector
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

        F = (alpha(rho) * inner(u, v) * dx + mu * inner(grad(u), grad(v)) * dx +
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

        def eval(self, _u, _v):
            _v.apply('insert')
            A_u = _v.copy()
            self.A.mult(_u, A_u)
            return _v.inner(A_u)

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
            if iteration >= 0:
                control_file << self.rho
                state_file << self.state

    class VolConstraint(ROL.Constraint):

        def __init__(self, inner_product):
            ROL.Constraint.__init__(self)
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
    x = FeVector(x.vector(), dot_product)
    g = Function(A)
    g = FeVector(g.vector(), dot_product)
    d = interpolate(Expression("1 + x[0] * (1-x[0])*x[1] * (1-x[1])", degree=1), A)
    d = FeVector(d.vector(), dot_product)
    x.checkVector(d, g)

    jd = Function(A)
    jd = FeVector(jd.vector(), dot_product)

    lower = interpolate(Constant(0.0), A)
    lower = FeVector(lower.vector(), dot_product)
    upper = interpolate(Constant(1.0), A)
    upper = FeVector(upper.vector(), dot_product)

    # Instantiate Objective class for poisson problem
    obj = ObjR(dot_product)
    # obj.checkGradient(x, d, 4, 2)
    volConstr = VolConstraint(dot_product)
    volConstr.checkApplyJacobian(x, d, jd, 3, 1)
    volConstr.checkAdjointConsistencyJacobian(v, d, x)

    set_log_level(30)
    
    paramsDict = {
            'General': {
                'Secant': { 'Type': 'Limited-Memory BFGS', 'Maximum Storage': 25 } },
                'Step': {
                    'Type': 'Augmented Lagrangian',
                    'Line Search': {
                        'Descent Method': {
                          'Type': 'Quasi-Newton Step'}},
                    'Augmented Lagrangian': {
                        'Initial Penalty Parameter'               : 1.e2,
                        'Penalty Parameter Growth Factor'         : 2,
                        'Minimum Penalty Parameter Reciprocal'    : 0.1,
                        'Initial Optimality Tolerance'            : 1.0,
                        'Optimality Tolerance Update Exponent'    : 1.0,
                        'Optimality Tolerance Decrease Exponent'  : 1.0,
                        'Initial Feasibility Tolerance'           : 1.0,
                        'Feasibility Tolerance Update Exponent'   : 0.1,
                        'Feasibility Tolerance Decrease Exponent' : 0.9,
                        'Print Intermediate Optimization History' : True,
                        'Subproblem Step Type'                    : 'Line Search',
                        'Subproblem Iteration Limit'              : 10
                      }},
            'Status Test': {
                'Gradient Tolerance': 1e-15, 'Relative Gradient Tolerance': 1e-10,
                'Step Tolerance': 1e-16, 'Relative Step Tolerance': 1e-10,
                'Iteration Limit': 7}
            }
    params = ROL.ParameterList(paramsDict, "Parameters")
    bound_constraint = ROL.Bounds(lower, upper, 1.0)

    optimProblem = ROL.OptimizationProblem(obj, x, bnd=bound_constraint, econ=volConstr, emul=l)
    solver = ROL.OptimizationSolver(optimProblem, params)
    solver.solve()
