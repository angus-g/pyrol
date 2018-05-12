import ROL


class HS4_Obj(ROL.Objective):
    def __init__(self):
        super().__init__()

    def value(self, x, tol):
        return ((x[0]+1)**3)/3. + x[1]

    def gradient(self, g, x, tol):
        g[0] = (x[0]+1)**2
        g[1] = 1


def HS4_Bnd(lower, upper):
    lower[0] = 1
    lower[1] = 0
    upper[0] = 1000
    upper[1] = 1000

def HS4_initial_guess(x):
    x[0] = 1.125
    x[1] = 0.125

def HS4_minimum(x):
    def near(a, b): return abs(a-b) < 1e-6
    return near(x[0], 1.0) and near(x[1], 0.)


def run_HS4(Vec, params_dict):
    obj = HS4_Obj()
    x = Vec(2)
    HS4_initial_guess(x)
    lower = Vec(2)
    upper = Vec(2)
    HS4_Bnd(lower, upper)
    bnd = ROL.Bounds(lower, upper, 1.0)
    params = ROL.ParameterList(params_dict, "Parameters")
    problem = ROL.OptimizationProblem(obj, x, bnd=bnd)
    solver = ROL.OptimizationSolver(problem, params)
    solver.solve()
    print(x[0], x[1])
    assert HS4_minimum(x)
