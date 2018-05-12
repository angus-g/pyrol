import ROL


class HS28_Obj(ROL.Objective):
    def __init__(self):
        super().__init__()

    def value(self, x, tol):
        return (x[0]+x[1])**2 + (x[1]+x[2])**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] + x[1])
        g[1] = 2 * (x[0] + x[1]) + 2 * (x[1] + x[2])
        g[2] = 2 * (x[1] + x[2])

    def hessVec(self, hv, v, x, tol):
        # 2     2   0
        # 2     4   2
        # 0     2   2
        hv[0] = 2 * v[0] + 2 * v[1] + 0*v[2]
        hv[1] = 2 * v[0] + 4 * v[1] + 2*v[2]
        hv[2] = 0 * v[0] + 2 * v[1] + 2*v[2]

class HS28_Econ(ROL.Constraint):
    def __init__(self):
        super().__init__()

    def value(self, cvec, x, tol):
        cvec[0] = x[0] + 2.0 * x[1] + 3.0 * x[2] - 1.0

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = v[0] + 2 * v[1] + 3 * v[2]
  
    def applyAdjointJacobian(self, jv, v, x, tol):
        jv[0] = v[0]
        jv[1] = 2 * v[0]
        jv[2] = 3 * v[0]

    def applyAdjointHessian(self, ahuv, u, v, x, tol):
        ahuv.scale(0.)


def HS28_initial_guess(x):
    x[0] = -4.
    x[1] = 1.
    x[2] = 1.

def HS28_minimum(x):
    def near(a, b): return abs(a-b) < 1e-6
    return near(x[0], 0.5) and near(x[1], -0.5) and near(x[2], 0.5)


def run_HS28(Vec, params_dict):
    obj = HS28_Obj()
    x = Vec(3)
    HS28_initial_guess(x)
    # obj.checkGradient(x)
    # obj.checkHessVec(x)
    HS28_initial_guess(x)
    l = Vec(1)
    l[0] = 0.0
    con = HS28_Econ()
    params = ROL.ParameterList(params_dict, "Parameters")
    problem = ROL.OptimizationProblem(obj, x, econs=[con], emuls=[l])
    solver = ROL.OptimizationSolver(problem, params)
    solver.solve()
    print(x[0], x[1], x[2])
    assert HS28_minimum(x)


if __name__ == "__main__":
    run_HS28(ROL.StdVector, {})
