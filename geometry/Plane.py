import numpy as np
from numpy.linalg import norm
import math

from geometry.Line import Line


class Plane:
    def __init__(self, normal, tolerance, p=None, plane=None):
        if plane is None:
            n = norm(normal)
            if n < 1e-10:
                raise Exception("Norm is zero!")
            self.w = normal / n
            self.tolerance = tolerance
            self.origin_offset = -np.dot(p, self.w) if p is not None else 0
            self.origin = -self.origin_offset * self.w
            self.u = self.orthogonal(self.w)
            self.v = np.cross(self.w, self.u)
        else:
            self.origin_offset = plane.origin_offset
            self.origin = plane.origin
            self.u = plane.u
            self.v = plane.v
            self.w = plane.w
            self.tolerance = plane.tolerance

    def orthogonal(self, w):
        threshold = 0.6 * norm(w)
        if threshold == 0:
            raise Exception("Norm is zero!")
        x = w[0]
        y = w[1]
        z = w[2]
        if abs(x) <= threshold:
            inverse = 1 / math.sqrt(y ** 2 + z ** 2)
            return np.array([0, inverse * z, -inverse * y])
        elif abs(y) <= threshold:
            inverse = 1 / math.sqrt(x ** 2 + z ** 2)
            return np.array([-inverse * z, 0, inverse * x])
        inverse = 1 / math.sqrt(x ** 2 + y ** 2)
        return np.array([inverse * y, -inverse * x, 0])

    def get_normal(self):
        return self.w

    def get_origin(self):
        return self.origin

    def project(self, p):
        p2 = np.array([np.dot(p, self.u), np.dot(p, self.v)])
        return self.u * p2[0] + self.v * p2[1] - self.origin_offset * self.w

    def get_point_at(self, in_plane, offset):
        return self.u * in_plane[0] + self.v * in_plane[1] - \
               self.origin_offset * self.w

    def intersection(self, l=None, other=None):
        if l is not None:
            dir = l.get_direction()
            dot = np.dot(self.w, dir)
            if dot < 1e-10:
                return None
            p = l.point_at(0.0)
            k = - (self.origin_offset + np.dot(self.w, p)) / dot
            return p + k * dir
        else:
            dir = np.cross(self.w, other.w)
            if norm(dir) < self.tolerance:
                return None
            p = self.intersection_3_planes(self, other, Plane(dir,
                                                            self.tolerance))
            return Line(p1=p, p2=p + dir, tolerance=self.tolerance)

    def intersection_3_planes(self, p1, p2, p3):
        a1 = p1.w[0]
        b1 = p1.w[1]
        c1 = p1.w[2]
        d1 = p1.origin_offset

        a2 = p2.w[0]
        b2 = p2.w[1]
        c2 = p2.w[2]
        d2 = p2.origin_offset

        a3 = p3.w[0]
        b3 = p3.w[1]
        c3 = p3.w[2]
        d3 = p3.origin_offset

        a23 = b2 * c3 - b3 * c2
        b23 = c2 * a3 - c3 * a2
        c23 = a2 * b3 - a3 * b2

        determinant = a1 * a23 + b1 * b23 + c1 * c23
        if abs(determinant) < 1e-10:
            return None

        r = 1.0 / determinant

        return np.array([
        (-a23 * d1 - (c1 * b3 - c3 * b1) * d2 - (c2 * b1 - c1 * b2) * d3) * r,
        (-b23 * d1 - (c3 * a1 - c1 * a3) * d2 - (c1 * a2 - c2 * a1) * d3) * r,
        (-c23 * d1 - (b1 * a3 - b3 * a1) * d2 - (b2 * a1 - b1 * a2) * d3) * r
        ])

    def contains(self, p):
        return abs(self.get_offset(point=p)) < self.tolerance

    def get_offset(self, point=None, plane=None):
        if plane is None:
            return np.dot(point, self.w) + self.origin_offset
        else:
            return self.origin_offset + (-plane.origin_offset if
                    self.same_orientation_as(plane) else plane.origin_offset)

    def same_orientation_as(self, other):
        return np.dot(self.w, other.w) > 0.0