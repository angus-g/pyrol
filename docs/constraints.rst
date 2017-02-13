Constraints
===========

ROL supports three different types of constraints.

Bound Constraints
-----------------

Creating a bound constraint is a simple as 

.. code-block:: python

    x_lo = ROL.StdVector(2)
    x_lo[0] = -1.0
    x_lo[1] = -1.0
    x_up = ROL.StdVector(2)
    x_up[0] = +1.0
    x_up[1] = +1.0
    scale = 1.0
    bnd = ROL.BoundConstraint(x_lo,x_up, scale)
    

Note that in order to use this with your own vector class, the class needs to implement the `__getitem__` method.

Equality Constraints
--------------------

In order to write your own equality constraint, inherit from the class `ROL.EqualityConstraint`

Inequality Constraints
----------------------

Not yet implemented in PyROL.
