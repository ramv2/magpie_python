from collections import OrderedDict
from numpy.linalg import norm
from Vassal.AtomImage import AtomImage
import numpy as np

class VoronoiCell:
    def __init__(self, atom, radical, faces=None):
        self.atom = atom
        self.faces = faces
        self.radical = radical
        self.volume = float("inf")

    def get_atom(self):
        return self.atom

    def get_faces(self):
        return self.faces

    def n_faces(self):
        return len(self.faces)

    def get_vertices(self):
        return set([face.get_vertices() for face in self.faces])

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
                    new_shell[AtomImage(neighbor.get_atom(), old_image)] = None

            # Remove all images corresponding to the shell inside this one.
            for shell in previous_shells:
                new_shell.pop(shell)

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

        # Now that all the paths are gathered, output only thr last step and
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
        if self.volume == float("inf"):
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

            if not set(self.faces) <= set(face.get_neighboring_faces()):
                return False

        return True