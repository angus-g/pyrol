import ROL

class FiredrakeVector(ROL.Vector):
    def __init__(self, vec, inner):
        ROL.Vector.__init__(self)
        self.vec = vec
        self.inner = inner

    def plus(self, x):
        self.vec += x.vec

    def scale(self, alpha):
        self.vec *= alpha

    def __getitem__(self, i):
        return self.vec[i]

    def __setitem__(self, i, v):
        self.vec[i] = v

    def reduce(self, r, r0):
        res = r0
        tempx = self.vec.get_local()
        for i in range(len(tempx)):
            res = r(tempx[i], res)
        self.vec.set_local(tempx)
        return res

    def applyBinary(self, f, y):
        tempx = self.vec.get_local()
        tempy = y.vec.get_local()
        for i in range(len(tempx)):
            tempx[i] = f(tempx[i], tempy[i])
        self.vec.set_local(tempx)

    def applyUnary(self, f):
        res = r0
        tempx = self.vec.get_local()
        for i in range(len(tempx)):
            tempx[i] = f(tempx[i])
        self.vec.set_local(tempx)

    def dot(self, xx):
        if self.inner is not None:
            return self.inner.eval(self.vec, xx.vec)
        else:
            return self.vec.inner(xx.vec)

    def clone(self):
        res = FiredrakeVector(self.vec.copy(), self.inner)
        return res

    def dimension(self):
        return len(self.vec)

    def basis(self, i):
        res = FiredrakeVector(self.vec.copy(), self.inner)
        res.scale(0.0)
        res[i] = 1.0
        return res
