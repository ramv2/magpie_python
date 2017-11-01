from numpy.linalg import norm
import numpy as np

class VoronoiVertex:
    def __init__(self, inside_atom=None, position=None, edge1=None, edge2=None):
        if position is not None and not isinstance(position, np.ndarray):
            raise ValueError("Position should be a numpy array!.")

        in_atom = inside_atom
        pos = position

        self.previous_edge = edge1
        self.next_edge = edge2

        if in_atom is None and pos is None:
            in_atom = edge1.get_edge_face().get_inside_atom()
            pos = edge1.get_line().intersection(edge2.get_line())

        self.position = pos
        self.distance = norm(self.position-np.array(
            in_atom.get_position_cartesian()))

        if inside_atom is None and position is None:
            if edge1.is_ccw(edge2=edge2):
                self.previous_edge = edge1
                self.next_edge = edge2
            else:
                self.previous_edge = edge2
                self.next_edge = edge1

    @classmethod
    def get_centroid(self, points):
        center = np.zeros(3, dtype=float)
        for p in points:
            center += p.get_position()

        center /= len(points)
        return center

    def distance_from(self, vertex):
        return norm(self.position - vertex.position)

    def get_distance_from_center(self):
        return self.distance

    def get_position(self):
        return self.position

    def __str__(self):
        return str(self.position)

    def __hash__(self):
        return 3 + hash(self.position)

    def __eq__(self, other):
        if isinstance(other, VoronoiVertex):
            return other.previous_edge.__eq__(self.previous_edge) and \
                   other.next_edge.__eq__(self.next_edge)
        return False

    def __cmp__(self, other):
        if self.previous_edge.__eq__(other.previous_edge):
            return self.next_edge.__cmp__(other.next_edge)
        else:
            return self.previous_edge.__cmp__(other.previous_edge)

    def get_next_edge(self):
        return self.next_edge

    def get_previous_edge(self):
        return self.previous_edge

    def is_on_edge(self, e):
        return self.previous_edge.__eq__(e) or self.next_edge.__eq(e)