Objectives
==========

To write your own objective class, simply inherit from `ROL.Objective`.

Let's say we want to minimize the function :math:`f(x,y) = (x-1)^2 + y^2`.

.. code-block:: python

    class MyObj(ROL.Objective):
        def __init__(self):
            ROL.Objective.__init__(self)

        def value(self, x, tol):
            return (x[0] - 1)**2 + x[1]**2

        def gradient(self, g, x, tol):
            g[0] = 2 * (x[0] - 1)
            g[1] = 2 * x[1]

If we omit the definition of the `gradient` function, then ROL will use a finite difference approximation instead.
