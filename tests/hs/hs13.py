import ROL


class HS13_Obj(ROL.Objective):
    def __init__(self):
        super().__init__()

    def value(self, x, tol):
        return (x[0] - 2.)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2. * (x[0]-2.)
        g[1] = 2. * x[1]
    
    def hessVec(self, hv, v, x, tol):
        hv[0] = 2. * v[0]
        hv[1] = 2. * v[1]


class HS13_Icon(ROL.Constraint):
    def __init__(self):
        super().__init__()

    def value(self, cvec, x, tol):
        cvec[0] = ((1.-x[0])**3)-x[1]

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = -3. * (1.-x[0])**2 * v[0] - v[1]
  
    def applyAdjointJacobian(self, jv, v, x, tol):
        jv[0] = -3. * (1.-x[0])**2 * v[0]
        jv[1] = -v[0]

    def applyAdjointHessian(self, ahuv, u, v, x, tol):
        # 6 * (1-x[0])  0
        # 0             0
        ahuv.scale(0.)
        ahuv[0] = u[0] * 6 * (1-x[0]) * v[0]

def HS13_Ibnds(lower):
    lower[0] = 0.0

def HS13_Bnd(lower):
    lower[0] = 0
    lower[1] = 0

def HS13_initial_guess(x):
    x[0] = -2.
    x[1] = -2.

def HS13_minimum(x):
    def near(a, b): return abs(a-b) < 1e-5
    return near(x[0], 1.0) and near(x[1], 0.)

def run_HS13(Vec, params_dict):
    obj = HS13_Obj()
    x = Vec(2)
    d = Vec(2)
    HS13_initial_guess(x)
    l = Vec(1)
    l[0] = 1.0
    icon = HS13_Icon()
    ilower = Vec(1)
    HS13_Ibnds(ilower)
    ibnd = ROL.Bounds(ilower, isLower=True)
    lower = Vec(2)
    HS13_Bnd(lower)
    bnd = ROL.Bounds(lower, isLower=True)
    params = ROL.ParameterList(params_dict, "Parameters")
    problem = ROL.OptimizationProblem(obj, x, bnd=bnd, icons=[icon], imuls=[l], ibnds=[ibnd])
    solver = ROL.OptimizationSolver(problem, params)
    solver.solve()
    print(x[0], x[1])
    assert HS13_minimum(x)

if __name__ == "__main__":
    run_HS13(ROL.StdVector, my_dict)
