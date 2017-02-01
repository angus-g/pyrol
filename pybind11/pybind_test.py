import ROL
vec = ROL.StdVector(2)
# for v in vec:
#     print v

# print vec.norm()
vec[0] = 0.0
vec[1] = 0.0
# for v in vec:
#     print v
# for v in vec:
#     print v
# print vec.norm()

class MyObj(ROL.Objective):
    def __init__(self):
        ROL.Objective.__init__(self)

    def value(self, x, tol):
        return (x[0] - 1)**2 + x[1]**2
    
    def gradient(self, g, x, tol):
        g[0] = 2 * (x[0] - 1)
        g[1] = 2 * x[1]

obj = MyObj()
print(obj.value(vec, 0))
g = ROL.StdVector(2)
g[0] = 0.0
g[1] = 0.0
obj.gradient(g, vec, 0.001)

print(vec[0])
print(vec[1])
algo = ROL.Algorithm("Line Search", {"test":"test"})
algo.run(vec, obj)
print(vec[0])
print(vec[1])

# print(g[0])
# print(g[1])
