import ROL
import numpy as np

class NumpyVector(ROL.Vector):
    def __init__(self, size):
        ROL.Vector.__init__(self)
        self.data = np.zeros(size, dtype=np.float64)

    def plus(self, xx):
        self.data += xx.data

    def scale(self, alpha):
        self.data *= alpha

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, v):
        self.data[i] = v

    def dot(self, xx):
        return np.inner(self.data, xx.data)

    def dimension(self):
        return len(self.data)

    def basis(self, i):
        res = NumpyVector(self.dimension())
        res[i] = 1.0
        return res

    def clone(self):
        res = NumpyVector(self.dimension())
        # res.data = np.copy(self.data)
        return res

    def applyBinary(self, f, x):
        for i in range(len(self.data)):
            self.data[i] = f(self.data[i], x.data[i])

    def applyUnary(self, f):
        for i in range(len(self.data)):
            self.data[i] = f(self.data[i])

    def reduce(self, r, r0):
        res = r0
        for i in range(len(self.data)):
            res = r(self.data[i], res)
        return res
