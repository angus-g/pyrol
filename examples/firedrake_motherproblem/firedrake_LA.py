import ROL


class FiredrakeLA(ROL.CustomLA):
    def __init__(self, vec, inner):
        ROL.CustomLA.__init__(self)
        self.vec = vec
        self.inner = inner

    def plus(self, xx):
        self.vec += xx.vec

    def scale(self, alpha):
        self.vec *= alpha

    def dot(self, xx):
        if self.inner is not None:
            return self.inner.eval(self.vec, xx.vec)
        else:
            return self.vec.inner(xx.vec)

    def clone(self):
        res = FiredrakeLA(self.vec.copy(), self.inner)
        return res
