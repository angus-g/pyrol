Poisson Motherproblem
=====================

We want to solve the motherproblem of PDE Constrained Optimization as described at `dolfin-adjoint <http://www.dolfin-adjoint.org/en/latest/documentation/poisson-mother/poisson-mother.html/>`_.
We can use either dolfin or firedrake to do this.

Begin by including either of the two finite element libraries and the corresponding vector class from ROL ::

    from firedrake import UnitSquareMesh, FunctionSpace, Function, \
        Expression, TrialFunction, TestFunction, inner, grad, DirichletBC, \
        dx, Constant, solve, assemble, as_backend_type, File, sin, \
        SpatialCoordinate
    from ROL.firedrake_vector import FiredrakeVector as FeVector
    backend = "firedrake"
    #  from ROL.dolfin_vector import DolfinVector as FeVector
    #  from dolfin import UnitSquareMesh, FunctionSpace, Function, \
    #      Expression, TrialFunction, TestFunction, inner, grad, DirichletBC, \
    #      dx, Constant, solve, assemble, as_backend_type, File, CellFunction, refine
    #  backend = "dolfin"

Import ROL ::

    import ROL

Define a mesh and in the case of dolfin, optionally perform some random mesh refinement in order to show of the mesh indenpendent behaviour ::

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

        k = 3
        for i in range(k):
           mesh = randomly_refine(mesh)
    else:
        k = 0

We can either treat the optimization as an optimization in R^d, essentially forgetting about the fact that we optimize in a function space or we can use the proper Riesz map and inner product.
This is needed for mesh indenpendent behaviour. ::

    use_correct_riesz_and_inner = True

Set an output directory ::

    outdir = "output_riesz_%s_refinements_%i/" % (use_correct_riesz_and_inner, k)

Define the function space, regularity parameter and target function ::

    V = FunctionSpace(mesh, "Lagrange", 1)  # space for state variable
    M = FunctionSpace(mesh, "DG", 0)  # space for control variable
    beta = 1e-4
    yd = Function(V)
    x, y = SpatialCoordinate(mesh)
    pi = 3.141592653589793
    yd.interpolate((1.0/(2*pi*pi)) * sin(pi*x) * sin(pi*y))
    # uncomment this for a 'more difficult' target distribution
    # yd.interpolate((x <= 0.5)*(y <= 0.5))


Define function that solves for the state given a control ::

    def solve_state(u):
        y = TrialFunction(V)
        v = TestFunction(V)
        a = inner(grad(y), grad(v)) * dx
        L = u * v * dx
        bc = DirichletBC(V, Constant(0.0), "on_boundary")
        y = Function(V)
        solve(a == L, y, bc)
        return y

In order to obtain the gradient of the reduced functional we need the solution of the adjoint equation ::

    def solve_adjoint(u, y):
        lam = TrialFunction(V)
        v = TestFunction(V)
        a = inner(grad(lam), grad(v)) * dx
        L = -(y-yd) * v * dx
        bc = DirichletBC(V, Constant(0.0), "on_boundary")
        lam = Function(V)
        solve(a == L, lam, bc)
        return lam

Define the proper inner product based on the L^2 inner product on the function space ::

    class L2Inner(object):

        def __init__(self):
            self.A = assemble(TrialFunction(M)*TestFunction(M)*dx)
            self.Ap = as_backend_type(self.A).mat()
            bcs = [DirichletBC(M, Constant(0.0), "on_boundary")]
            self.Ap_b = assemble(TrialFunction(M)*TestFunction(M)*dx, bcs=bcs)

        def eval(self, _u, _v):
            upet = as_backend_type(_u).vec()
            vpet = as_backend_type(_v).vec()
            A_u = self.Ap.createVecLeft()
            self.Ap.mult(upet, A_u)
            return vpet.dot(A_u)

        def riesz_map(self, derivative):
            if backend == "firedrake":
                rhs = Function(M, val=derivative.dat)
                res = Function(M)
                solve(self.Ap_b, res, rhs)
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

Define output files ::

    state_file = File(outdir + "state.pvd")
    control_file = File(outdir + "control.pvd")

Define the objective class, inheriting from ROL.Objective ::

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
            if iteration >= 0:
              if backend == "firedrake":
                  control_file.write(self.u)
                  state_file.write(self.y)
              else:
                  control_file << self.u
                  state_file << self.y

Set some basic parameters for the optimization. We want to use L-BFGS for the optimization ::

    paramsDict = {
        "Step": {
            "Line Search": {
                "Descent Method": {
                    "Type": "Quasi-Newton Method"
                }
            },
            "Type": "Line Search",
        },
        "Status Test": {
            "Gradient Tolerance": 1e-10,
            "Iteration Limit": 20
        }
    }
    params = ROL.ParameterList(paramsDict, "Parameters")

Create the inner product ::

    if use_correct_riesz_and_inner:
        inner_product = L2Inner()
    else:
        inner_product = None

Create the objective ::

    obj = Objective(inner_product)

Create vectors for the optimization and perform a linear algebra check::

    u = Function(M)
    opt = FeVector(u.vector(), inner_product)
    d = Function(M)
    d.interpolate(sin(x*pi)*sin(y*pi))
    d = FeVector(d.vector(), inner_product)
    # if backend == "firedrake":
    #     obj.checkGradient(opt, d, 3, 1)

Create the upper and lower bound constraints ::

    xlo = Function(M)
    xlo.interpolate(Constant(0.0))
    x_lo = FeVector(xlo.vector(), inner_product)
    xup = Function(M)
    xup.interpolate(Constant(0.9))
    x_up = FeVector(xup.vector(), inner_product)
    bnd = ROL.Bounds(x_lo, x_up, 1.0)

Run the optimization ::

    algo = ROL.Algorithm("Line Search", params)
    algo.run(opt, obj, bnd)
    if backend == "firedrake":
        File("res.pvd").write(Function(M, opt.vec))
    else:
        File("res.pvd") << Function(M, opt.vec)
