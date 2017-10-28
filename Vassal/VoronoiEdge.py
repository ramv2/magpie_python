import numpy as np
from numpy.linalg import norm
from sympy import Line

from Vassal.VoronoiVertex import VoronoiVertex


class VoronoiEdge:
    def __init__(self, edge_face=None, intersecting_face=None):

        self.edge_face = edge_face
        self.intersecting_face = intersecting_face
        self.line = None
        self.direction = None
        self.beginning = -float("inf")
        self.next_edge = None
        self.end = float("inf")
        self.previous_edge = None

        # Compute the line.
        self.line = edge_face.get_plane().intersection(
            intersecting_face.get_plane())[0]
        if not self.line:
            raise Exception("Planes are parallel.")

        # Ensure vector is CCW w.r.t edge face.
        cut_direction = -intersecting_face.get_normal()
        self.direction = np.array(self.line.direction.unit.evalf(), dtype=float)
        if not self.is_ccw(vec1=edge_face.get_normal(),
                           vec2=self.direction, vec3=cut_direction):
            self.line = Line(self.line.p1, -self.line.p2)

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

    @classmethod
    def get_abscissa(self, point, line):
        p1 = np.array(line.p1.evalf(), dtype=float)
        p2 = np.array(line.p2.evalf(), dtype=float)
        zero = p1 - np.dot(p1, p2 - p1) * (p2 - p1)/norm(p2 - p1)
        direction = np.array(line.direction.unit.evalf(), dtype=float)
        return np.dot(point - zero, direction)

    @classmethod
    def compute_intersection(self, edge1, edge2, just_join=False):
        # Determine the point at which the edges intersect.
        point = np.array(edge1.line.intersection(edge2.line)[0].evalf(),
                         dtype=float)
        if point is None or point.size == 0:
            if just_join:
                raise Exception("Edges do not intersect.")
            else:
                return False

        # Determine the relationship between edges (using their directions).
        is_forward = edge1.is_ccw(edge2=edge2)

        # Using the direction, check whether intersection is within bounds of
        #  each edge.
        edge1_terminus = self.get_abscissa(point, edge1.get_line())
        edge2_terminus = self.get_abscissa(point, edge2.get_line())

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
        h = 43 * h + id(self.edge_face)
        h = 43 * h + id(self.intersecting_face)
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

    def find_next_edge(self, candidates):
        # Locate the ccw edges.
        ccw_edges = []
        for edge in candidates:
            if self.is_ccw(edge2=edge):
                ccw_edges.append(edge)

        # Check if any were found.
        if len(ccw_edges) == 0:
            return None
        elif len(ccw_edges) == 1:
            return ccw_edges[0]

        # Find the closest edge(s).
        closest_edges = []
        min_dist = float("inf")
        min_point = np.array([np.inf]*3)
        for edge in ccw_edges:
            other_line = edge.get_line()
            # If this line contains the minimum point, add it to list.
            if not np.isinf(min_point).any() and other_line.contains(min_point):
                closest_edges.append(edge)
                continue

            intersection = self.line.intersection(other_line)
            if not intersection:
                # Line is anti parallel.
                continue
            intersection = np.array(intersection[0].evalf(), dtype=float)

            # See if it is the closest.
            x = self.get_abscissa(intersection, self.line)
            if x < min_dist:
                closest_edges = []
                min_dist = x
                min_point = intersection
                closest_edges.append(edge)

        # If only one edge, return answer.
        if len(closest_edges) == 1:
            return closest_edges[0]

        # Otherwise, find the edge with the largest angle.
        max_angle = 0
        choice = None
        for edge in closest_edges:
            angle = self.line.angle_between(edge.get_line()).evalf()
            if angle > max_angle:
                choice = edge
                max_angle = angle

        return choice

    def get_start_vertex(self):
        return VoronoiVertex(edge1=self, edge2=self.previous_edge)

    def get_end_vertex(self):
        return VoronoiVertex(edge1=self, edge2=self.next_edge)

    def generate_pair(self):
        try:
            return VoronoiEdge(self.get_intersecting_face(),
                               self.get_edge_face())
        except Exception:
            raise Exception("Shouldn't be possible.")