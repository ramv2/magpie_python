import numpy as np
import math
from decimal import Decimal
from numpy.linalg import norm
from Vassal.VoronoiEdge import VoronoiEdge
from Vassal.VoronoiVertex import VoronoiVertex
from geometry.Plane import Plane


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
        self.tol = 1e-8

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
        h = 19 + h * id(self.inside_atom)
        h = 19 + h * id(self.outside_atom)
        return h

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
            inside_pos = self.inside_atom.get_position_cartesian()
            outside_pos = self.outside_atom.get_position()
            self.face_plane = Plane(p=self.face_center, normal=
                                    outside_pos-inside_pos, tolerance=1e-8)
        return self.face_plane

    def get_normal(self):
        if self.face_normal is None:
            self.face_normal = self.get_plane().get_normal()
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
                next_vertex = self.vertices[(i + 1) % l]
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
        plane = self.get_plane()
        if plane.contains(point):
            return 0
        offset = plane.get_offset(point=point)
        if offset > 0:
            return 1
        else:
            return -1

    def __str__(self):
        output = str(self.inside_atom.get_id()+"->"+self.outside_atom.__str__())
        return output

    def assemble_face_from_faces(self, faces):
        # Generate a local copy of this list.
        other_faces = list(faces)
        if self in other_faces:
            other_faces.remove(self)

        # Find all edges.
        available_edges = []
        for other_face in other_faces:
            edge = None
            try:
                edge = VoronoiEdge(self, other_face)
            except Exception:
                continue

            # Store the data.
            if edge is None:
                continue
            else:
                available_edges.append(edge)

        # Error check.
        if len(available_edges) < 3:
            raise Exception("Not enough edges.")

        # for i,e in enumerate(available_edges):
        #     l = e.get_line()
        #     print i
        #     print l.zero
        #     print l.direction
        #     print l.distance(p=self.face_center)
            # print i,e.intersecting_face.outside_atom
        # print "Face center: "+str(self.face_center)
        cur_edge = available_edges[0]
        for i, edge in enumerate(available_edges):
            flag = cur_edge.is_ccw(edge2=edge)
            print i, edge.intersecting_face.outside_atom, flag

        available_edges.sort(cmp=self.edge_comp)

        return self.assemble_face_from_edges(available_edges)

    def edge_comp(self, a, b):
        d1 = a.get_line().distance(p=self.face_center)
        d2 = b.get_line().distance(p=self.face_center)

        # l1 = a.get_line()
        # print "A values: " + a.intersecting_face.outside_atom.__str__()
        # print l1.zero
        # print l1.direction
        # print l1.tolerance
        # print d1

        # l2 = b.get_line()
        # print "B values: " + b.intersecting_face.outside_atom.__str__()
        # print l2.zero
        # print l2.direction
        # print l2.tolerance
        # print d2
        #
        # print
        # print
        if d1 < d2:
            return -1
        elif d1 > d2:
            return 1
        else:
            return 0

    def assemble_face_from_edges(self, available_edges):
        # Clear cached face area.
        self.face_area = np.nan

        # The closest edge is certainly on the face.
        face_edges = []
        face_edges.append(available_edges[0])
        cur_edge = face_edges[0]

        # Now, loop through all the edges and compute their intersection with
        # all other edges on the face. This is accomplished by finding the
        # next edge and proceeding until we reach the beginning or none is
        # found.
        while True:
            next_edge = cur_edge.find_next_edge(available_edges)
            if next_edge == face_edges[0]:
                break
            elif next_edge is None:
                raise Exception("Face is not closed.")
            elif len(face_edges) > len(available_edges):
                # This happens if the first vertex is not a valid edge.
                # find the first edge and continues from a different starting
                # point. To address this, easy way is to remove the problem
                # vertex and try again.
                new_avail = list(available_edges)
                new_avail.remove(new_avail[0])
                if len(new_avail) < 3:
                    raise Exception("Face failed to build.")
                return self.assemble_face_from_edges(new_avail)
            face_edges.append(next_edge)
            cur_edge = next_edge

        # Compute intersections and add edges.
        self.edges = []
        l_fe = len(face_edges)
        for i in range(l_fe):
            VoronoiEdge.compute_intersection(face_edges[i], face_edges[(i +
                                    1) % l_fe], just_join=True)
            self.edges.append(face_edges[i])

        # Store vertices.
        self.compute_vertices()

        return self.is_closed()


    def compute_vertices(self):
        self.vertices = []
        for e in self.edges:
            try:
                self.vertices.append(e.get_start_vertex())
            except Exception:
                raise Exception("Vertex computation error.")

        # Remove area.
        self.get_area()

    def is_closed(self):
        if len(self.edges) < 3:
            return False

        for i in range(len(self.vertices)-1):
            if not (self.edges[i].get_next_edge() in self.edges and
                        self.edges[i].get_previous_edge() in self.edges):
                return False
        return True

    def compute_intersection(self, new_face):
        # Clear cached area.
        self.face_area = np.nan

        # See if any vertex is outside this new face or if any vertices are
        # on this face.
        outside_vertex = None
        n_c = 0
        for v in self.vertices:
            rel_pos = new_face.position_relative_to_face(v.get_position())
            if rel_pos > 0:
                outside_vertex = v
                break
            elif rel_pos == 0:
                n_c += 1

        # If it doesn't intersect, return "None".
        if outside_vertex is None:
            # Check if two vertices are on the new face. If so, they must
            # share an edge because these faces are convex.
            if n_c < 2:
                return None

            # If so, find what edge they correspond to.
            new_edge = None
            edge_to_replace = None
            try:
                new_edge = VoronoiEdge(self, new_face)
            except Exception:
                raise Exception("New face does not intersect this face.")

            try:
                for e in self.edges:
                    if new_face.position_relative_to_face(e.get_start_vertex(

                    ).get_position()) == 0 and \
                                    new_face.position_relative_to_face(
                                        e.get_end_vertex().get_position()) ==\
                                    0:
                        edge_to_replace = e
                        break
            except Exception:
                raise Exception("Face was not closed to begin with.")

            # Compute the new edges.
            if edge_to_replace is None:
                raise Exception("Edge on the new face was not found.")

            try:
                VoronoiEdge.compute_intersection(new_edge,
                    edge_to_replace.get_previous_edge(), just_join=True)
                VoronoiEdge.compute_intersection(new_edge,
                    edge_to_replace.get_next_edge(), just_join=True)
            except Exception:
                raise Exception("New edge and old are not parallel.")

            # Replace edge in the lineup.
            self.edges[self.edges.index(edge_to_replace)] = new_edge

            # Recompute vertices.
            self.compute_vertices()

            return new_edge.generate_pair()

        # Get this vertex ID number.
        outside_vertex_id = self.vertices.index(outside_vertex)

        # Traverse forward until a vertex inside is found.
        last_vertex = outside_vertex_id
        l_v = len(self.vertices)
        while new_face.position_relative_to_face(self.vertices[
                                    last_vertex].get_position()) >= 0:
            last_vertex = (last_vertex + 1) %  l_v
            # If we go all the way around the loop, all vertices are outside
            # this face. Clear everything and return "None" because no edge
            # is formed.
            if outside_vertex_id == last_vertex:
                self.edges = []
                self.vertices = []
                return None

        # Traverse backwards until you find the first vertex that is inside
        # the intersecting plane.
        first_vertex = outside_vertex_id
        while new_face.position_relative_to_face(self.vertices[
                                    first_vertex].get_position()) >= 0:
            first_vertex -= 1
            if first_vertex < 0:
                first_vertex += l_v

        # Create the new edge.
        new_edge = None
        try:
            new_edge = VoronoiEdge(self, new_face)
        except Exception:
            raise Exception("New face doesn't intersect this face.")

        start_edge = self.vertices[first_vertex].get_next_edge()
        end_edge = self.vertices[last_vertex].get_previous_edge()
        # Eliminate edges that are outside the face.
        cur_edge = start_edge.get_next_edge()
        while cur_edge != end_edge:
            if cur_edge not in self.edges:
                raise Exception("Edge not found.")
            else:
                self.edges.remove(cur_edge)

            cur_edge = cur_edge.get_next_edge()

        # Join start and edges with this new edge.
        try:
            VoronoiEdge.compute_intersection(start_edge, new_edge,
                                             just_join=True)
            VoronoiEdge.compute_intersection(end_edge, new_edge, just_join=True)
        except Exception:
            raise Exception("New edge doesn't intersect current ones.")

        # Add it to the edges list.
        idx = self.edges.index(start_edge)
        self.edges.insert(1 + idx, new_edge)
        # Recompute vertices.
        self.compute_vertices()
        return new_edge.generate_pair()

    def reset_edges(self, edges):
        try:
            l_e = len(edges)
            for i in range(l_e):
                VoronoiEdge.compute_intersection(edges[i], edges[ (i + 1) %
                                            l_e], just_join=True)
            self.edges = edges
            self.compute_vertices()
        except Exception:
            raise RuntimeError("Was the stored state in error?")


    def get_cut_length(self, new_face):
        # Find a single vertex that is outside the face.
        outside_vertex = None
        for v in self.vertices:
            if new_face.position_relative_to_face(v.get_position()) > 0:
                outside_vertex = v
                break

        # If no vertex is outside, cut length == 0.
        if outside_vertex is None:
            return 0

        # Get the vertex ID number.
        outside_vertex_id = self.vertices.index(outside_vertex)

        # Traverse forward until a vertex inside is found.
        last_vertex = outside_vertex_id
        l_v = len(self.vertices)
        while new_face.position_relative_to_face(self.vertices[
                                                     last_vertex].get_position()) >= 0:
            last_vertex = (last_vertex + 1) % l_v
            # If we go all the way around the loop, all vertices are outside
            # this face. Clear everything and return "None" because no edge
            # is formed.
            if outside_vertex_id == last_vertex:
                self.edges = []
                self.vertices = []
                return float("inf")

        # Traverse backwards until you find the first vertex that is inside
        # the intersecting plane.
        first_vertex = outside_vertex_id
        while new_face.position_relative_to_face(self.vertices[
                                                     first_vertex].get_position()) >= 0:
            first_vertex -= 1
            if first_vertex < 0:
                first_vertex += l_v

        # Create the new edge.
        new_edge = None
        try:
            new_edge = VoronoiEdge(self, new_face)
        except Exception:
            raise Exception("New face doesn't intersect this face.")

        start_edge = self.vertices[first_vertex].get_next_edge()
        end_edge = self.vertices[last_vertex].get_previous_edge()

        p1 = start_edge.get_line().intersection(new_edge.get_line())
        p2 = end_edge.get_line().intersection(new_edge.get_line())
        return norm(p1 - p2)

    def is_contacted_by(self, other_face):
        for v in self.vertices:
            if other_face.position_relative_to_face(v.get_position()) >= 0:
                return True
        return False

    def is_completely_outside(self, other_face):
        for v in self.vertices:
            if other_face.position_relative_to_face(v.get_position()) < 0:
                return False
        return True

    def contains_vertex(self, v):
        return v in self.vertices