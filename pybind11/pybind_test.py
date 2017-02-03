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

# x = ROL.StdVector(2)
# x[0] = -1.0
# x[1] = 2.0

# print x[0] 
# print x[1] 
# algo.run(x, obj)
# print x[0] 
# print x[1] 


# _x = ROL.EigenVector(2)
# _x[0] = -1.0
# _x[1] = 2.0

# print _x[0] 
# print _x[1] 
# algo = ROL.Algorithm("Line Search", {"test":"test"})
# algo.run(_x, obj)
# print _x[0] 
# print _x[1] 

# import numpy as np
# a = np.asarray([0, 1, 2, 3])
# evec = ROL.EigenVector(a)
# print evec[0]
# print evec[1]
# print evec[2]
# print evec[3]
# a[3] = 5.0
# print evec[3]
# evec[3]=3.1
# print a[3]

# import ROL
import numpy as np
def iplus(x, y):
    x.data[0] += y.data[0]
    x.data[1] += y.data[1]

def iscale(x, alpha):
    x.data[0] *= alpha

class NPBasedLA(ROL.CustomLA):
    def __init__(self):
        ROL.CustomLA.__init__(self)
        # super(ROL.CustomLA, self).__init__(self)
        self.data = np.asarray([0, 0])
        self.clones = []
        print "construct"

    def iplus(self, xx):
        print "here plus"

    def iscale(self, alpha):
        print "here scale"

    def iclone(self):
        print "here clone"
        res = NPBasedLA()
        self.clones.append(res)
        # res.iscale(1.0)
        return res

x = NPBasedLA()
y = NPBasedLA()
z = NPBasedLA()
# x.scale(3.0)
# x.plus(y)
# # xx = x.clone()
# # xx.scale(1.0)
# yy = x.iclone()
# yy.iscale(1.0)
x.checkClone()
x.checkVector(y, z)
# # algo.run(x, obj)
