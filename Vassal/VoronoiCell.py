from collections import OrderedDict
from numpy.linalg import norm
from Vassal.AtomImage import AtomImage
import numpy as np
from Vassal.VoronoiEdge import VoronoiEdge
from Vassal.VoronoiFace import VoronoiFace

class VoronoiCell:
    def __init__(self, atom, radical, faces=None):
        self.atom = atom
        self.faces = faces
        self.radical = radical
        self.volume = np.nan

    def get_atom(self):
        return self.atom

    def get_faces(self):
        return self.faces

    def n_faces(self):
        return len(self.faces)

    def get_vertices(self):
        return list(set([face.get_vertices() for face in self.faces]))

    def get_neighbor_types(self):
        output = OrderedDict()
        for face in self.faces:
            atom_id = face.get_outside_atom().get_atom_id()
            atom_type = self.atom.get_cell().get_atom(atom_id).get_type()
            if atom_type in output:
                output[atom_type] += 1
            else:
                output[atom_type] = 1

        return output

    def get_neighbors(self):
        return [face.get_outside_atom() for face in self.faces]

    def get_neighbor_distances(self):
        return [face.get_neighbor_distance() for face in self.faces]

    def get_neighbor_shells(self, cells, index):
        if index < 0:
            # Special case: Bad input.
            raise ValueError("Shell index must be >= 0.")
        elif index == 0:
            # Special case: 0th cell (this atom).
            return [OrderedDict({AtomImage(self.atom, [0, 0, 0]):None})]
        else:
            # General case: i-th shell.
            # Get the previous shell.
            previous_shells = self.get_neighbor_shells(cells, index - 1)

            # For each atom in the previous shell, get its neighbors.
            new_shell = OrderedDict()
            for atom in previous_shells[index - 1]:
                # Get its nearest neighbor according to the diagram.
                neighbors = cells[atom.get_atom_id()].get_neighbors()
                # For each of these neighbors, adjust the image coordinates.
                this_image = atom.get_supercell()
                for neighbor in neighbors:
                    # Get the coordinates relative to the central atom.
                    old_image = neighbor.get_supercell() + this_image

                    # Create a new image, and append it to output.
                    new_image = AtomImage(neighbor.get_atom(), old_image)
                    if new_image not in new_shell:
                        new_shell[new_image] = None

            # Remove all images corresponding to the shell inside this one.
            for shell in previous_shells:
                for image in shell:
                    if image in new_shell:
                        new_shell.pop(image)

            # Append new shell to output.
            previous_shells.append(new_shell)

            return previous_shells

    def get_extended_faces(self, cells, index):
        if index == 0:
            # Special case: Index == 0.
            return list(self.faces)

        # Get neighbor list.
        neighbor_list = self.get_neighbor_shells(cells, index)

        # Last shell are atoms on the outside.
        outside_images = neighbor_list[index]

        # Get atoms to consider.
        all_images = OrderedDict({[[image for images in neighbor_list for
                                    image in images]]:None})

        # Get the faces of all outside atoms that do not correspond to an
        # atom on the inside of the polyhedron.
        output = []
        for image in outside_images:
            # Get coordinates of this image w.r.t. central atom.
            image_coord = image.get_supercell()
            # Get its Voronoi cell.
            cell = cells[image.get_atom_id()]
            # Iterate over its faces.
            for face in cell.get_faces():
                # Get the image on the outside of this face.
                relative_image = face.get_outside_atom()
                # Get the coordinates of this image relative to the central
                # atom.
                outside_coord = relative_image.get_supercell() + image_coord

                # Get the actual image.
                actual_image = AtomImage(relative_image.get_atom(),
                                         outside_coord)

                # If this image is not inside of the polyhedron, add the face
                #  to the output.
                if actual_image not in all_images:
                    output.append(face)
        return output

    def get_neighbors_by_walks(self, cells, shell):
        # Gather all possible paths.
        # Initialize with first step.
        paths = []
        # Add in starting point.
        paths.append(([AtomImage(self.atom, [0, 0, 0])], 1.0))
        # Increment until desired shell.
        for step in range(shell):
            # Get the new list.
            new_paths = []
            # For each current path.
            for path in paths:
                # Get last atom in current path.
                last_step = path[0][len(path[0]) - 1]
                # Get each possible step.
                new_steps = OrderedDict()
                # Surface area of possible steps.
                surface_area = 0.0
                for face in cells[last_step.get_atom_id()].get_faces():
                    # Adjust image based on position of last atom.
                    new_supercell = face.get_outside_atom().get_supercell() +\
                                    last_step.get_supercell()

                    # Store the next step.
                    next_step = AtomImage(face.get_outside_atom().get_atom(),
                                          new_supercell)
                    area = face.get_area()
                    surface_area += area
                    new_steps[next_step] = area

                # Eliminate backtracking steps.
                for previous_step in path[0]:
                    if previous_step in new_steps:
                        surface_area -= new_steps.pop(previous_step)

                # Create new paths, making sure to update weights.
                for k,v in new_steps.iteritems():
                    # Increment path.
                    new_path = list(path[0])
                    new_path.append(k)

                    # Increment weight.
                    new_weight = path[1] * v / surface_area

                    # Add it to new_paths.
                    new_paths.append((new_path, new_weight))

            # Update paths.
            paths = new_paths

        # Now that all the paths are gathered, output only the last step and
        # weights of all paths that lead to that step.
        output = OrderedDict()
        for path in paths:
            # Get the last step.
            atom = path[0][len(path[0]) - 1]

            # Update map.
            if atom in output:
                output[atom] += path[1]
            else:
                output[atom] = path[1]

        return output

    def get_neighbor_shell(self, cells, index):
        return self.get_neighbor_shells(cells, index)[index]

    def get_num_shared_bonds(self, cell, direction, neighbors):
        n_shared = 0
        for face in cell.get_faces():
            other_cell_n_id = face.get_outside_atom().get_atom_id()
            other_cell_n_image = face.get_outside_atom().get_supercell() + \
                                 direction
            # Make that image.
            other_cell_n = AtomImage(self.atom.get_cell().get_atom(
                other_cell_n_id), other_cell_n_image)

            # Check against every neighbor of this cell.
            for n in neighbors:
                if other_cell_n.__eq__(n):
                    n_shared += 1

        return n_shared

    def get_coordination_shell_shape(self, cells):
        # Get IDs of every neighbor.
        neighbors = [face.get_outside_atom() for face in self.faces]

        # Get number of mutual neighbors for each neighbor.
        output = OrderedDict()
        for n in neighbors:
            n_shared = self.get_num_shared_bonds(cells[n.get_atom_id()],
                                                 n.get_supercell(), neighbors)
            if n_shared in output:
                output[n_shared] += 1
            else:
                output[n_shared] = 1

        return output

    def get_polyhedron_shape(self):
        output = OrderedDict()
        for face in self.faces:
            n_vertices = face.n_vertices()
            if n_vertices in output:
                output[n_vertices] += 1
            else:
                output[n_vertices] = 1

        return output

    def get_volume(self):
        if np.isnan(self.volume):
            self.volume = 0
            atom_center = self.atom.get_position_cartesian()
            for face in self.faces:
                area = face.get_area()
                from_center = face.get_centroid() - atom_center
                n = face.get_normal()
                n /= norm(n)
                h = np.dot(from_center, n)
                self.volume += area * h / 3.0

        return self.volume

    def get_surface_area(self):
        return sum([face.get_area() for face in self.faces])

    def get_min_max_vertex_distance(self):
        vertices = self.get_vertices()
        l_v = len(vertices)
        min_dist = float("inf")
        max_dist = -min_dist
        for i in range(l_v):
            for j in range(i+1, l_v):
                dist = vertices[i].distance_from(vertices[j])
                min_dist = min(dist, min_dist)
                max_dist = max(dist, max_dist)

        return min_dist, max_dist

    def geometry_is_valid(self):
        for face in self.faces:
            if not face.is_closed():
                return False

            for f in face.get_neighboring_faces():
                if f not in self.faces:
                    return False
        return True

    def compute_cell(self, image_finder, cutoff):
        cur_cutoff = cutoff
        n_attempts = 0
        while n_attempts < 4:
            n_attempts += 1
            image_finder.set_cutoff_distance(cur_cutoff)
            # Find all nearby images.
            images = [image[0] for image in
                      image_finder.get_all_neighbors_of_atom(
                          self.atom.get_id())]

            # Compute cell.
            try:
                self.compute_cell_helper(images)
            except Exception:
                cur_cutoff *= 1.5
                continue

            return

        raise Exception("Cell failed to compute.")

    def compute_cell_helper(self, images):
        # Clear cached volume.
        self.volume = np.nan

        # Get all possible faces.
        possible_faces = self.compute_faces(images)

        # Get the faces corresponding to the direct polyhedron.
        direct_faces = self.compute_direct_neighbors(possible_faces)

        # Construct direct polyhedron.
        for df in direct_faces:
            try:
                df.assemble_face_from_faces(direct_faces)
            except Exception:
                raise Exception("Direct polyhedron failed to construct.")

        self.faces = list(direct_faces)

        # Get the faces that might actually be direct faces.
        possible_indirect_faces = self.compute_possible_indirect_neighbors(
            possible_faces)

        # Use these faces to compute indirect neighbors.
        for face in possible_indirect_faces:
            self.compute_intersection(face)

    def compute_faces(self, images):
        # Generate faces.
        output = []
        for image in images:
            # Check if the image is this atom.
            if image.get_atom_id() == self.atom.get_id() and np.array_equal(
                    image.get_supercell(), [0, 0, 0]):
                # If so, skip.
                continue

            # Make the appropriate face.
            try:
                output.append(VoronoiFace(self.atom, image, self.radical))
            except Exception:
                raise RuntimeError("Failure to create face.")

        return output

    def compute_direct_neighbors(self, faces):
        # Sort distance from face to the center.
        faces.sort(cmp=self.compare_faces)

        # The closest face is on the direct polyhedron.
        direct_faces = []
        direct_faces.append(faces.pop(0))

        # Now loop through all faces to find those whose centers are not
        # outside any other face that is on the direct polyhedron.

        to_remove = []
        for face in faces:
            is_inside = True
            for df in direct_faces:
                d = df.position_relative_to_face(face.get_face_center())
                if d >= 0:
                    is_inside = False
                    break

            if is_inside:
                direct_faces.append(face)
                to_remove.append(face)
        for face in to_remove:
            faces.remove(face)
        return direct_faces

    def compare_faces(self, a, b):
        d1 = a.get_face_distance()
        d2 = b.get_face_distance()
        if d1 < d2:
            return -1
        elif d1 > d2:
            return 1
        else:
            return 0

    def compute_possible_indirect_neighbors(self, possible_faces):
        # Get list of faces that might be indirect neighbors.
        max_vertex_distance = max([v.get_distance_from_center() for face in
                                   self.faces for v in face.get_vertices()])

        possible_indirect_faces = []
        for face in possible_faces:
            if face.get_face_distance() < max_vertex_distance:
                possible_indirect_faces.append(face)
            else:
                # possible faces is sorted by get_direct_faces.
                break

        return possible_indirect_faces

    def compute_intersection(self, new_face):
        # Clear cached result.
        self.volume = np.nan

        # Determine whether any faces intersect this new new_face.
        no_intersection = True
        for cur_face in self.faces:
            for v in cur_face.get_vertices():
                d = new_face.position_relative_to_face(v.get_position())
                if d > 0:
                    no_intersection = False
                    break

        if no_intersection:
            return False

        # If it does, store the old face information.
        previous_state = OrderedDict({face: face.get_edges() for face in
                                     self.faces})

        # Attempt to perform intersections.
        try:
            new_edges = []
            to_remove = []
            for c_face in self.faces:
                new_edge = c_face.compute_intersection(new_face)
                if c_face.n_edges() < 3:
                    to_remove.append(c_face)
                else:
                    if new_edge is not None:
                        new_edges.append(new_edge)
            for f in to_remove:
                self.faces.remove(f)
            # Check if we have enough edges.
            if len(new_edges) < 3:
                raise Exception("Not enough edges were formed.")

            # Assemble new face and add it to list of faces.
            new_face.assemble_face_from_edges(new_edges)
            self.faces.append(new_face)

            # Check if geometry is valid.
            if not self.geometry_is_valid():
                raise Exception("Geometry is invalid.")
            return True
        except Exception:
            # Restore previous state.
            self.faces = []
            for face in previous_state:
                face.reset_edges(previous_state[face])
                self.faces.append(face)
            return False

    def remove_face(self, to_remove):
        if to_remove not in self.faces:
            raise RuntimeError("No such face exists.")

        self.faces.remove(to_remove)
        self.volume = np.nan

        # Find all faces that are currently in contact with this face.
        contacting_faces = to_remove.get_neighboring_faces()
        for face_to_build in contacting_faces:
            cur_edges = face_to_build.get_edges()

            # Compute edges corresponding to the "contacting faces".
            for face in contacting_faces:
                if face.__eq__(face_to_build):
                    continue
                new_edge = VoronoiEdge(face_to_build, face)
                cur_edges.append(new_edge)

            # Remove the edge corresponding to the face being removed.
            for e in cur_edges:
                if e.get_intersecting_face().__eq__(to_remove):
                    cur_edges.remove(e)
            face_to_build.assemble_face_from_edges(cur_edges)