#!/usr/bin/env py.test

import ROL

# Create a quadratic objective function around (1.0, 0.0)
class Objective(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2

    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]
obj = Objective()

def test_std_vector_check():
    x = ROL.StdVector(2)
    y = ROL.StdVector(2)
    z = ROL.StdVector(2)

    x[0] = 1.5
    x[1] = 0.5
    y[0] = 1.2
    y[1] = 0.2

    u = x.checkVector(y, z)
    assert sum(u) < 1e-12


def test_std_vector_run():

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
    params = ROL.ParameterList(params_dict, "Parameters")
    algo = ROL.Algorithm("Line Search", params)
    x = ROL.StdVector(2)
    x[0] = -1.0
    x[1] = 2.0

    algo.run(x, obj)
    # Check answer to 8 decimal places
    assert round(x[0] - 1.0, 8) == 0.0
    assert round(x[1], 8) == 0.0
