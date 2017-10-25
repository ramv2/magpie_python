from sympy import Plane, Point3D
from numpy.linalg import norm
import numpy as np
from Vassal.VoronoiEdge import VoronoiEdge
from Vassal.VoronoiVertex import VoronoiVertex

class VoronoiFace:
    def __init__(self, inside_atom, outside_atom, radical):
        self.inside_atom = inside_atom
        self.outside_atom = outside_atom
        self.face_plane = None
        self.face_normal = None
        self.edges = []
        self.vertices = []
        self.face_distance = None
        self.face_center = None
        self.face_area = np.nan

        inside_pos = self.inside_atom.get_position_cartesian()
        outside_pos = self.outside_atom.get_position()
        diff = inside_pos - outside_pos
        atom_dist = norm(diff)
        if radical:
            self.face_distance = self.get_plane_distance(
                self.inside_atom.get_radius(), self.outside_atom.get_atom(
                ).get_radius(), atom_dist)
        else:
            self.face_distance = atom_dist / 2

        self.face_center = inside_pos - diff * self.face_distance / atom_dist

    def get_plane_distance(self, r1, r2, d):
        return (r1**2 - r2**2 + d**2)/ 2 / d

    def __eq__(self, other):
        if isinstance(other, VoronoiFace):
            return other.inside_atom.get_id() == self.inside_atom.get_id() \
                   and other.outside_atom.__eq__(self.outside_atom)
        return False

    def __hash__(self):
        h = 7
        h = 19 + h * hash(self.inside_atom)
        h = 19 + h * hash(self.outside_atom)

    def __cmp__(self, other):
        if self.inside_atom.get_id() == other.inside_atom.get_id():
            return self.outside_atom.__cmp__(other.outside_atom)
        else:
            return self.inside_atom.get_id() - other.inside_atom.get_id()

    def get_inside_atom(self):
        return self.inside_atom

    def get_face_center(self):
        return self.face_center

    def get_outside_atom(self):
        return self.outside_atom

    def get_plane(self):
        if self.face_plane is None:
            inside_pos = self.inside_atom.get_position_cartesian
            outside_pos = self.outside_atom.get_position()
            self.face_plane = Plane(self.face_center,
                                    normal_vector=outside_pos-inside_pos)
        return self.face_plane

    def get_normal(self):
        if self.face_plane is None:
            self.face_normal = np.array(self.get_plane().normal_vector)
        return self.face_normal

    def n_edges(self):
        return len(self.edges)

    def get_area(self):
        # If needed compute area.
        if np.isnan(self.face_area):
            # Get centroid of face.
            centroid = self.get_centroid()
            # Loop over all edges.
            area = 0.0
            l = len(self.vertices)
            for i in range(l):
                this_vertex = self.vertices[i]
                next_vertex = self.vertices[i + 1 % l]
                a = this_vertex.get_position() - centroid
                b = next_vertex.get_position() - centroid
                area += norm(np.cross(a, b))

            self.face_area = area / 2

        return self.face_area

    def get_centroid(self):
        return VoronoiVertex.get_centroid(self.vertices)

    def get_face_distance(self):
        return self.face_distance

    def get_neighbor_distance(self):
        return norm(self.inside_atom.get_position_cartesian() -
                    self.outside_atom.get_position())

    def get_vertices(self):
        return list(self.vertices)

    def n_vertices(self):
        return len(self.vertices)

    def get_edges(self):
        return list(self.edges)

    def get_neighboring_faces(self):
        output = []
        for e in self.edges:
            f = e.get_intersecting_face()
            if f not in output:
                output.append(f)

        return output

    def get_common_vertices(self, other_face):
        output = set(self.vertices)
        return output.intersection(other_face.vertices)

    def position_relative_to_face(self, point):
        if self.face_plane.is_coplanar(Point3D(point)):
            return 0

        offset = np.dot(point, self.get_normal())
        if offset > 0:
            return 1
        else:
            return -1

    def __str__(self):
        output = str(self.inside_atom.get_id()+"->"+self.outside_atom.__str__())
        return output

    def assemble_from_vertices(self, vertices):
        # Store the vertices.
        self.vertices = []
        self.vertices = vertices

        # Create the edges.
        self.edges = []
        l = len(self.vertices)
        for i in range(l):
            edge = VoronoiEdge(start_vertex=self.vertices[i],
                               end_vertex=self.vertices[i + 1 % l],
                               edge_face=self)
            self.edges.append(edge)

    def is_closed(self):
        return len(self.vertices) >= 3