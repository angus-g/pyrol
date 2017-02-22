import ROL
from dolfin import as_backend_type

class dolfinLA(ROL.CustomLA):
    def __init__(self, vec, inner=None):
        ROL.CustomLA.__init__(self)
        self.vec = vec
        self.inner = inner

    def plus(self, x):
        # PETSc doesn't like adding vector to themselves,
        # so we have to add this check
        xpet = as_backend_type(x.vec).vec()
        thispet = as_backend_type(self.vec).vec()
        if xpet == thispet:
            self.vec *= 2
        else:
            self.vec += x.vec

    def scale(self, alpha):
        self.vec *= alpha

    def __getitem__(self, i):
        return self.vec[i][0]

    def __setitem__(self, i, v):
        self.vec[i] = v

    def dot(self, xx):
        if self.inner is not None:
            return self.inner.eval(self.vec, xx.vec)
        else:
            return self.vec.inner(xx.vec)

    def clone(self):
        res = dolfinLA(self.vec.copy(), self.inner)
        return res

    def dimension(self):
        return len(self.vec)

    def basis(self, i):
        res = dolfinLA(self.vec.copy(), self.inner)
        res.scale(0.0)
        res[i] = 1.0
        return res
