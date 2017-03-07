import ROL

class FiredrakeLA(ROL.CustomLA):
    def __init__(self, vec, inner):
        ROL.CustomLA.__init__(self)
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

    def dot(self, xx):
        if self.inner is not None:
            return self.inner.eval(self.vec, xx.vec)
        else:
            return self.vec.inner(xx.vec)

    def clone(self):
        res = FiredrakeLA(self.vec.copy(), self.inner)
        return res

    def dimension(self):
        return len(self.vec)

    def basis(self, i):
        res = FiredrakeLA(self.vec.copy(), self.inner)
        res.scale(0.0)
        res[i] = 1.0
        return res
