from numpy.linalg import norm
import numpy as np

class VoronoiVertex:
    def __init__(self, position, inside_atom=None, faces=None):
        if not isinstance(position, np.ndarray):
            raise ValueError("Position should be a numpy array!.")

        if inside_atom is None and faces is None:
            raise ValueError("Needs either inside_atom or list of faces.")

        self.position = position
        in_atom = inside_atom
        if in_atom is None:
            in_atom = faces[0].get_inside_atom()

        self.distance = norm(self.position-np.array(
            in_atom.get_position_cartesian()))

        self.faces = set(faces)

    def distance_from(self, vertex):
        return norm(self.position - vertex.position)

    def get_distance_from_center(self):
        return self.distance

    def get_position(self):
        return self.position

    def __hash__(self):
        return self.faces.__hash__()

    @classmethod
    def get_centroid(self, points):
        center = np.zeros(3, dtype=float)
        for p in points:
            center += p.get_position()

        center /= len(points)
        return center

    def __eq__(self, other):
        if isinstance(other, VoronoiVertex):
            return other.faces.__eq__(self.faces)
        return False

    def get_faces(self):
        return self.faces

    def __cmp__(self, other):
        # If face sizes are different.
        if len(other.faces) != len(self.faces):
            return len(self.faces) - len(other.faces)

        # Make sorted list of the faces.
        my_faces = list(sorted(self.faces))
        your_faces = list(sorted(other.faces))

        # Compare each entry.
        for i in range(len(my_faces)):
            comp = my_faces[i].__cmp__(your_faces[i])
            if comp != 0:
                return comp

        return 0