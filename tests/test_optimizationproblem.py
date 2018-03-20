import ROL
params_dict = {
    'General': {
        'Secant': {
            'Type': 'Limited-Memory BFGS',
            'Maximum Storage': 5
        }
    },
    'Step': {
        'Type': 'Line Search',
        'Line Search': {
            'Descent Method': {
                'Type': 'Quasi-Newton Method'
            },
            'Curvature Condition': {
                'Type': 'Strong Wolfe Conditions'
            }
        }
    },
    'Status Test': {
        'Gradient Tolerance': 1e-15,
        'Relative Gradient Tolerance': 1e-10,
        'Step Tolerance': 1e-16,
        'Relative Step Tolerance': 1e-10,
        'Iteration Limit': 10
    }
}
class Objective(ROL.Objective):
    def __init__(self):
        super().__init__()

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]

obj = Objective()
x = ROL.StdVector(2)
x[0] = -1.0
x[1] = 2.0

lower = ROL.StdVector(2)
lower[0] = -10
lower[1] = -10

upper = ROL.StdVector(2)

upper[0] = 0.5
upper[1] = 0.5
bnd = ROL.Bounds(lower, upper, 1.0)

params = ROL.ParameterList(params_dict, "Parameters")
problem = ROL.OptimizationProblem(obj, x)
solver = ROL.OptimizationSolver(problem, params)
solver.solve()

# algo = ROL.Algorithm("Line Search", params)
# algo.run(x, obj)

print(x[0])
print(x[1])
problem = ROL.OptimizationProblem(obj, x, bnd=bnd)
solver = ROL.OptimizationSolver(problem, params)
solver.solve()
print(x[0])
print(x[1])
