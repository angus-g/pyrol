import ROL


class HS29_Obj(ROL.Objective):
    def __init__(self):
        super().__init__()

    def value(self, x, tol):
        return -x[0] * x[1] * x[2]

    def gradient(self, g, x, tol):
        g[0] = -x[1] * x[2]
        g[1] = -x[0] * x[2]
        g[2] = -x[0] * x[1]

    def hessVec(self, hv, v, x, tol):
        hv[0] = 0 * v[0] - x[2] * v[1] - x[1]*v[2]
        hv[1] = -x[2] * v[0] + 0 * v[1] -x[0]*v[2]
        hv[2] = -x[1] * v[0]  - x[0] * v[1] + 0*v[2]


class HS29_Icon(ROL.Constraint):
    def __init__(self):
        super().__init__()

    def value(self, cvec, x, tol):
        cvec[0] = -x[0] ** 2 - 2.0 * x[1] ** 2 - 4.0 * x[2] ** 2 + 48.0

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = -2.0 * x[0] * v[0] - 4.0 * x[1] * v[1] - 8.0 * x[2] * v[2]
  
    def applyAdjointJacobian(self, jv, v, x, tol):
        jv[0] = -2.0 * x[0] * v[0]
        jv[1] = -4.0 * x[1] * v[0]
        jv[2] = -8.0 * x[2] * v[0]

    def applyAdjointHessian(self, ahuv, u, v, x, tol):
        ahuv.scale(0.)
        ahuv[0] =-u[0] * (2*v[0])
        ahuv[1] =-u[0] * (4*v[1])
        ahuv[2] =-u[0] * (8*v[2])


def HS29_Ibnds(lower):
    lower[0] = 0.0


def HS29_initial_guess(x):
    x[0] = 1.
    x[1] = 1.
    x[2] = 1.


def HS29_minimum(x):
    a = 4
    b = 2*(2**0.5)
    c = 2
    def near(a, b): return abs(a-b) < 1e-5
    var1 = near(x[0], a) and near(x[1], b) and near(x[2], c)
    var2 = near(x[0], a) and near(x[1], -b) and near(x[2], -c)
    var3 = near(x[0], -a) and near(x[1], b) and near(x[2], -c)
    var4 = near(x[0], -a) and near(x[1], -b) and near(x[2], c)
    return var1 or var2 or var3 or var4


def run_HS29(Vec, params_dict):
    params = ROL.ParameterList(params_dict, "Parameters")
    obj = HS29_Obj()
    x = Vec(3); x[0] = 1.; x[1] = 2.; x[2] = .1
    d = Vec(3); d[0] = 1;  d[1] = -1; d[2] = 1.
    v = Vec(3); v[0] = 1;  v[1] = -1; v[2] = 1.
    # obj.checkGradient(x)
    # obj.checkHessVec(x, d, 4, 1)

    HS29_initial_guess(x)
    l = Vec(1)
    l[0] = 0.0
    con = HS29_Icon()
    jv = Vec(1); jv[0] = 1.
    # con.checkApplyJacobian(x, d, jv, 4, 1)
    # con.checkAdjointConsistencyJacobian(jv, d, x)
    # con.checkApplyAdjointHessian(x, jv, d, v, 5, 1)
    ilower = Vec(1)
    HS29_Ibnds(ilower)
    ibnd = ROL.Bounds(ilower, isLower=True)
    problem = ROL.OptimizationProblem(obj, x, icons=[con], imuls=[l], ibnds=[ibnd])
    solver = ROL.OptimizationSolver(problem, params)
    solver.solve()
    print(x[0], x[1], x[2])
    assert HS29_minimum(x)
