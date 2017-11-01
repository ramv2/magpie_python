import numpy as np
from numpy.linalg import norm
import sys

class Line:
    def __init__(self, p1=None, p2=None, tolerance=None, l=None):

        if l is None:
            delta = p2 - p1
            norm1 = norm(p2 - p1)
            norm2 = norm1 ** 2
            if norm2 == 0.0:
                raise Exception("Norm is zero!")
            self.direction = delta / norm1
            self.zero = p1 - np.dot(p1, delta) / norm2
            self.tolerance = 1e-10 if tolerance is None else tolerance
        else:
            self.tolerance = l.tolerance
            self.direction = l.direction
            self.zero = l.zero

    def set_tolerance(self, tol):
        self.tolerance = tol

    def get_tolerance(self):
        return self.tolerance

    def revert(self):
        reverted = Line(l=self)
        reverted.direction *= -1.0
        return reverted

    def get_direction(self):
        return self.direction

    def get_origin(self):
        return self.zero

    def get_abscissa(self, p):
        return np.dot(p - self.zero, self.direction)

    def point_at(self, abscissa):
        return self.zero + abscissa * self.direction

    def contains(self, p):
        return self.distance(p) < self.tolerance

    def distance(self, p=None, l=None):
        if l is None:
            d = p - self.zero
            n = d - np.dot(d, self.direction) * self.direction
            return norm(n)
        else:
            normal = np.cross(self.direction, l.direction)
            n = norm(normal)
            if n < sys.float_info.min:
                # Lines are parallel.
                return self.distance(p=l.zero)
            offset = np.dot(l.zero - self.zero, normal) / n
            return np.abs(offset)

    def closest_point(self, l):
        cos = np.dot(self.direction, l.direction)
        n = 1 - cos ** 2
        if n < sys.float_info.epsilon:
            return self.zero

        d0 = l.zero - self.zero
        a = np.dot(d0, self.direction)
        b = np.dot(d0, l.direction)
        return self.zero + self.direction * ( a - b * cos) / n

    def intersection(self, l):
        closest = self.closest_point(l)
        return closest if self.contains(closest) else None