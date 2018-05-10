#!/usr/bin/env py.test
import ROL
import numpy as np
from ROL.numpy_vector import NumpyVector
import pytest

class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]

class MyObj2(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return -(x[0]+x[1] + x[0]*x[1])

    def gradient(self, g, x, tol):
        g[0] = -(1 + x[1]);
        g[1] = -(1 + x[0]);


class EqConstraint(ROL.Constraint):
    def __init__(self):
        ROL.Constraint.__init__(self)

    def value(self, c, x, tol):
        c[0] = x[0]*x[0]+x[1]*x[1]-1

    def applyJacobian(self, jv, v, x, tol):
        jv[0] = 2*x[0]*v[0]+2*x[1]*v[1]

    def applyAdjointJacobian(self, ajv, v, x, tol):
        ajv[0] = v[0] * 2 * x[0]
        ajv[1] = v[0] * 2 * x[1]


paramsDict = {
    "Step": {
        "Line Search": {
            "Descent Method": {
                "Type": "Quasi-Newton Method"
            }
        },
        "Type": "Interior Point",
        # "Type": "Line Search",
        "Interior Point": {
            "Initial Barrier Penalty": 0.01,
            "Barrier Penalty Reduction Factor": 0.15,
            "Minimum Barrier Penalty": 1e-8
        },
        "Composite Step": {
            "Output Level": 0
        }
    },
    "Status Test": {
        "Gradient Tolerance" :1e-12,
        "Step Tolerance" :1e-16,
        "Constraint Tolerance" :1e-12,
        "Iteration Limit" :10
    }
}


def test_checkVectorPassed():
    x = NumpyVector(2)
    x.data[0] = 0.5
    x.data[1] = 0.5
    y = NumpyVector(2)
    y.data[0] = 0.2
    y.data[1] = 0.3
    z = NumpyVector(2)
    z.data[0] = 0.1
    z.data[1] = 0.7

    u = x.checkVector(y, z)
    assert sum(u) < 1e-12

def run_U(algo):
    obj = MyObj()
    paramsDict["Step"]["Type"] = algo
    params = ROL.ParameterList(paramsDict, "Parameters")
    x = NumpyVector(2)
    optimProblem = ROL.OptimizationProblem(obj, x)
    solver = ROL.OptimizationSolver(optimProblem, params)
    solver.solve()
    print(x.data)
    assert round(x[0] - 1.0, 6) == 0.0
    assert round(x[1], 6) == 0.0

def test_U_LS():
    run_U("Line Search")

def test_U_TR():
    run_U("Trust Region")

def test_U_BU():
    run_U("Bundle")

def run_B(algo):
    obj = MyObj()
    paramsDict["Step"]["Type"] = algo
    params = ROL.ParameterList(paramsDict, "Parameters")
    x = NumpyVector(2)
    x_lo = NumpyVector(2)
    x_lo[0] = -1
    x_lo[1] = -1
    x_up = NumpyVector(2)
    x_up[0] = +0.7
    x_up[1] = +0.7
    bnd = ROL.Bounds(x_lo, x_up, 1.0)
    optimProblem = ROL.OptimizationProblem(obj, x, bnd=bnd)
    solver = ROL.OptimizationSolver(optimProblem, params)
    solver.solve()
    assert round(x[0] - 0.7, 6) == 0.0
    assert round(x[1], 6) == 0.0

def test_B_LineSearch():
    run_B("Line Search")

def test_B_PDAS():
    run_B("Primal Dual Active Set")

def test_B_TR():
    run_B("Trust Region")

def test_B_MY():
    run_B("Moreau-Yosida Penalty")

@pytest.mark.xfail(reason="There is some issue with the interior point solver we need to figure out")
def test_B_IP():
    run_B("Interior Point")


def test_EqualityConstraintSatisfiesChecks():
    con = EqConstraint()
    x = NumpyVector(2)
    x[0] = 0.5 * 0.5**2
    x[1] = 0.5 * 0.5**2
    v = NumpyVector(2)
    v[0] = 1.0
    v[1] = -0.5
    jv = NumpyVector(1)
    w = NumpyVector(1)
    w[0] = 1.0
    con.checkApplyJacobian(x, v, jv, 4, 1);
    con.checkAdjointConsistencyJacobian(w, v, x)

def run_E(algo):
    obj = MyObj2()
    paramsDict["Step"]["Type"] = algo
    params = ROL.ParameterList(paramsDict, "Parameters")
    x = NumpyVector(2)
    x[0] = 0.5 * 0.5**2
    x[1] = 0.5 * 0.5**2
    l = NumpyVector(1)
    con = EqConstraint()
    optimProblem = ROL.OptimizationProblem(obj, x, econ=con, emul=l)
    solver = ROL.OptimizationSolver(optimProblem, params)
    solver.solve()
    assert round(x[0] - 0.707106, 5) == 0.0
    assert round(x[1] - 0.707106, 5) == 0.0

def test_E_AL():
    run_E("Augmented Lagrangian")

def test_E_CS():
    run_E("Composite Step")

def createBounds():
    x_lo = NumpyVector(2)
    x_lo[0] = -1
    x_lo[1] = -1
    x_up = NumpyVector(2)
    x_up[0] = +0.7
    x_up[1] = +0.7
    bnd = ROL.Bounds(x_lo, x_up, 1.0)
    bnd.test()
    return bnd

def test_create_bounds_seperately():
    obj = MyObj()
    paramsDict["Step"]["Type"] = "Trust Region"
    params = ROL.ParameterList(paramsDict, "Parameters")
    x = NumpyVector(2)
    bnd = createBounds()
    bnd.test()
    optimProblem = ROL.OptimizationProblem(obj, x, bnd=bnd)
    solver = ROL.OptimizationSolver(optimProblem, params)
    solver.solve()
    assert round(x[0] - 0.7, 6) == 0.0
    assert round(x[1], 6) == 0.0

def get_problem():
    obj = MyObj()
    x = NumpyVector(2)
    x_lo = NumpyVector(2)
    x_lo[0] = -1
    x_lo[1] = -1
    x_up = NumpyVector(2)
    x_up[0] = +0.7
    x_up[1] = +0.7
    bnd = ROL.Bounds(x_lo, x_up, 1.0)
    optimProblem = ROL.OptimizationProblem(obj, x, bnd=bnd)
    return optimProblem

def test_create_problem_seperately():
    paramsDict["Step"]["Type"] = "Trust Region"
    params = ROL.ParameterList(paramsDict, "Parameters")
    problem = get_problem()
    solver = ROL.OptimizationSolver(problem, params)
    solver.solve()
