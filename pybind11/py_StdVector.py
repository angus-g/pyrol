import ROL
class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2
    
    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]
obj = MyObj()
algo = ROL.Algorithm("Line Search", {"test":"test"})
x = ROL.StdVector(2)
x[0] = -1.0
x[1] = 2.0

print x[0] 
print x[1] 
algo.run(x, obj)
print x[0] 
print x[1] 
