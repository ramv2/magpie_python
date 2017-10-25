import numpy as np
from numpy.linalg import norm
from sympy import Line


class VoronoiEdge:
    def __init__(self, intersecting_face=None, start_vertex=None,
                 end_vertex=None, edge_face=None):

        if edge_face is None:
            raise ValueError("Edge face should be provided.")

        self.edge_face = edge_face
        self.intersecting_face = intersecting_face
        self.line = None
        self.direction = None
        self.beginning = -float("inf")
        self.next_edge = None
        self.end = float("inf")
        self.previous_edge = None
        self.start_vertex = start_vertex
        self.end_vertex = end_vertex

        if intersecting_face is not None:
            # Compute the line.
            self.line = edge_face.get_plane().intersection(
                intersecting_face.get_plane())
            if not self.line:
                raise RuntimeError("Planes are parallel.")

            # Ensure vector is CCW w.r.t edge face.
            cut_direction = -intersecting_face.get_normal()
            self.direction = np.array(self.line.direction.unit.evalf())
            if not self.is_ccw(vec1=edge_face.get_normal(),
                               vec2=self.direction, vec3=cut_direction):
                self.line = Line(self.line.p1, -self.line.p2)
        else:
            if start_vertex is None or end_vertex is None:
                raise ValueError("Start vertex and end vertex should be "
                                 "provided.")

            # Determine the intersecting face.
            face_candidates = set(start_vertex.get_faces())
            face_candidates.intersection(end_vertex.get_faces())
            face_candidates.remove(edge_face)
            if len(face_candidates) != 1:
                raise RuntimeError("Found wrong number of shared faces: "
                                   ""+str(len(face_candidates)))
            self.intersecting_face = face_candidates[0]

    def is_ccw(self, vec1=None, vec2=None, vec3=None, edge2=None):
        v1 = vec1
        v2 = vec2
        v3 = vec3
        if edge2 is not None:
            v1 = self.edge_face.get_normal()
            v2 = self.direction
            v3 = edge2.direction

        ccw = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) - v1[1] * (v2[0] * v3[
            2] - v2[2] * v3[0]) + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0])
        return ccw > 0

    def get_abscissa(self, point, edge):
        p1 = np.array(edge.line.p1)
        p2 = np.array(edge.line.p2)
        zero = p1 - np.dot(p1, p2 - p1) * (p2 - p1)/norm(p2 - p1)
        return np.dot(point - zero, edge.line.direction)


    def compute_intersection(self, edge1, edge2, just_join=False):
        # Determine the point at which the edges intersect.
        point = np.array(edge1.line.intersection(edge2.line))
        if not point:
            if just_join:
                raise ValueError("Edges do not intersect.")
            else:
                return False

        # Determine the relationship between edges (using their directions).
        is_forward = edge1.is_ccw(edge2=edge2)

        # Using the direction, check whether intersection is within bounds of
        #  each edge.
        edge1_terminus = self.get_abscissa(point, edge1)
        edge2_terminus = self.get_abscissa(point, edge2)

        if not just_join:
            within_bounds = None
            if is_forward:
                within_bounds = edge1_terminus < edge1.end and edge2_terminus > \
                                                               edge2.beginning
            else:
                within_bounds = edge1_terminus > edge1.beginning and \
                                edge2_terminus < edge2.end

            if not within_bounds:
                return False

        if is_forward:
            edge1.end = edge1_terminus
            edge1.next_edge = edge2
            edge2.beginning = edge2_terminus
            edge2.previous_edge = edge1
        else:
            edge1.beginning = edge1_terminus
            edge1.previous_edge = edge2
            edge2.end = edge2_terminus
            edge2.next_edge = edge1

        if not just_join:
            return True

    def __eq__(self, other):
        if isinstance(other, VoronoiEdge):
            return other.edge_face.__eq__(self.edge_face) and \
                   other.intersecting_face.__eq__(self.intersecting_face)
        return False

    def __hash__(self):
        h = 5
        h = 43 * h + hash(self.edge_face)
        h = 43 * h + hash(self.intersecting_face)
        return h

    def __cmp__(self, other):
        if self.edge_face.__eq__(other.edge_face):
            return self.intersecting_face.__cmp__(other.intersecting_face)
        return self.edge_face.__cmp__(other.edge_face)

    def get_line(self):
        return self.line

    def get_edge_face(self):
        return self.edge_face

    def get_intersecting_face(self):
        return self.intersecting_face

    def get_next_edge(self):
        return self.next_edge

    def get_previous_edge(self):
        return self.previous_edge

    def get_length(self):
        return norm((self.beginning - self.end) * self.direction)

    def __str__(self):
        output = "("+self.edge_face.__str__()+"," \
                        ""+self.intersecting_face.__str__()+")"
        return output

    def get_end_vertex(self):
        return self.end_vertex

    def get_start_vertex(self):
        return self.start_vertex