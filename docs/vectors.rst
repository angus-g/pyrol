Vectors
=======

To write your own vector class, simply inherit from `ROL.CustomLA`.

PyROL provides a `ROL.StdVector` class, which is based on the `ROL::StdVector` class.
Using it is as simple as 

.. code-block:: python

    import ROL
    x = ROL.StdVector(2)
    x[0] = 1.0
    x[1] = 2.0
    print x.norm()

However, you are likely to want to implement your own vector class, for example based on the data storage of a finite element library. 
This class is then able to implement a custom inner product, for example an :math:`L^2` inner product.


*TODO*: document one of the CustomLA classes and somehow include them here.
