import ROL

class dolfinLA(ROL.Vector):
    def __init__(self, vec, inner=None):
        ROL.Vector.__init__(self)
        self.vec = vec
        self.inner = inner

    def plus(self, x):
        # PETSc doesn't like adding a vector to itself,
        # so we have to make a copy. Would be good to fix this
        # upstream in DOLFIN
        self.vec += x.vec.copy()

    def scale(self, alpha):
        self.vec *= alpha

    def __getitem__(self, i):
        return self.vec[i][0]

    def __setitem__(self, i, v):
        # self.vec[i] = v
        a = self.vec.get_local()
        a[i] = v
        self.vec.set_local(a)
        # Total hack - apply when on the last entry
        if (i == self.vec.local_size() - 1):
            self.vec.apply("insert")

    def dot(self, xx):
        if self.inner is not None:
            return self.inner.eval(self.vec, xx.vec)
        else:
            return self.vec.inner(xx.vec)

    def clone(self):
        res = dolfinLA(self.vec.copy(), self.inner)
        return res

    def dimension(self):
        return self.vec.local_size()#len(self.vec)

    def basis(self, i):
        res = self.clone()
        res.scale(0.0)
        res[i] = 1.0
        return res
