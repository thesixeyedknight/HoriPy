# /home/sarthak/hori2/hori_project/hori/spatial_utils.py
import math

def create_spatial_grid(atoms_dict, cell_size):
	"""
	Partitions the 3D space occupied by atoms into a grid of cubic cells
	and assigns each atom to its corresponding cell.

	Args:
		atoms_dict (dict): Dictionary mapping atom ID to Atom namedtuple/object.
						   Each Atom object must have 'id', 'x', 'y', 'z' attributes.
		cell_size (float): The side length of each cubic cell.

	Returns:
		tuple: (grid, min_coords, cell_size)
			grid (dict): Keys are (ix, iy, iz) cell indices, values are lists of atom IDs.
			min_coords (tuple): (min_x, min_y, min_z) of the molecular bounding box.
			cell_size (float): The cell size used.
	"""
	if not atoms_dict or cell_size <= 0:
		return {}, (0.0, 0.0, 0.0), cell_size

	# 1. Determine Bounding Box
	# Initialize with coordinates of the first atom, or appropriate sentinels
	first_atom = next(iter(atoms_dict.values()))
	min_x, max_x = first_atom.x, first_atom.x
	min_y, max_y = first_atom.y, first_atom.y
	min_z, max_z = first_atom.z, first_atom.z

	for atom in atoms_dict.values():
		if atom.x < min_x: min_x = atom.x
		if atom.x > max_x: max_x = atom.x
		if atom.y < min_y: min_y = atom.y
		if atom.y > max_y: max_y = atom.y
		if atom.z < min_z: min_z = atom.z
		if atom.z > max_z: max_z = atom.z
	
	min_coords = (min_x, min_y, min_z)

	# 2. Initialize Grid
	grid = {}

	# 3. Assign Atoms to Cells
	for atom_id, atom in atoms_dict.items():
		# Ensure coordinates are relative to the grid's min_coords origin
		# Add a small epsilon for atoms exactly on min_coords to ensure they get index 0
		# if cell_size is chosen such that an atom on max_coord could be exactly on a boundary.
		# However, simple floor usually works.
		ix = math.floor((atom.x - min_x) / cell_size)
		iy = math.floor((atom.y - min_y) / cell_size)
		iz = math.floor((atom.z - min_z) / cell_size)
		
		# Handle atoms that might be exactly on the max boundary, they should belong to the last cell
		# Example: if max_x = 10, min_x = 0, cell_size = 5. Atom at x=10. (10-0)/5 = 2. floor(2) = 2.
		# Max cell index would be floor((max_x - min_x) / cell_size).
		# If atom.x == max_x, ix will be floor((max_x - min_x)/cell_size).
		# This could be one cell beyond the last cell containing any part of an atom *before* max_x.
		# A common approach is to subtract a tiny epsilon if an atom is exactly on max_x,
		# or ensure cell indices are capped.
		# For simplicity here, math.floor will assign it. If max_x = 10.0, atom.x = 10.0,
		# (10.0 - 0.0) / 5.0 = 2.0. ix = 2. Cells are 0, 1. Index 2 is correct for the start of cell 2.
		# If atom.x = 9.9, (9.9-0)/5 = 1.98. ix = 1.

		cell_key = (ix, iy, iz)
		if cell_key not in grid:
			grid[cell_key] = []
		grid[cell_key].append(atom_id) # Store atom ID

	return grid, min_coords, cell_size

def get_atoms_in_cutoff_from_grid(target_atom_coords, target_atom_id, all_atoms_dict, grid_data, cutoff_radius):
	"""
	Efficiently finds all atoms within a cutoff radius of a target atom's coordinates,
	using a pre-computed spatial grid.

	Args:
		target_atom_coords (tuple): (x, y, z) of the central/target atom.
		target_atom_id: The ID of the target atom itself (to exclude from neighbors).
		all_atoms_dict (dict): The main dictionary mapping atom ID to Atom object.
		grid_data (dict): Contains {'grid': grid_map, 'min_coords': (min_x, min_y, min_z), 'cell_size': size}.
		cutoff_radius (float): The search radius.

	Returns:
		list: A list of Atom objects (neighbors) within the cutoff_radius.
	"""
	if not grid_data or not all_atoms_dict or cutoff_radius <= 0:
		return []

	grid_map = grid_data['grid']
	min_c = grid_data['min_coords'] # min_coords is a tuple (min_x, min_y, min_z)
	cs = grid_data['cell_size']     # cell_size is a float

	found_neighbors = []

	# 1. Determine Target Atom's Cell Index
	tcx = math.floor((target_atom_coords[0] - min_c[0]) / cs)
	tcy = math.floor((target_atom_coords[1] - min_c[1]) / cs)
	tcz = math.floor((target_atom_coords[2] - min_c[2]) / cs)

	# 2. Determine Search Range of Cells based on cutoff_radius
	# How many cells in each direction (+/-) we need to check from the target cell.
	# If cutoff_radius is 7A and cell_size is 5A, ceil(7/5) = 2.
	# This means we check cells -2, -1, 0, +1, +2 relative to the target cell's index.
	cell_search_span = math.ceil(cutoff_radius / cs)

	# 3. Iterate Through Candidate Cells
	for ix_offset in range(-cell_search_span, cell_search_span + 1):
		for iy_offset in range(-cell_search_span, cell_search_span + 1):
			for iz_offset in range(-cell_search_span, cell_search_span + 1):
				candidate_cell_idx = (tcx + ix_offset, tcy + iy_offset, tcz + iz_offset)

				# 4. Check Atoms in Candidate Cell
				if candidate_cell_idx in grid_map:
					atom_ids_in_cell = grid_map[candidate_cell_idx]
					for atom_id in atom_ids_in_cell:
						if atom_id == target_atom_id: # Don't include the target atom itself
							continue

						neighbor_atom_obj = all_atoms_dict.get(atom_id)
						if not neighbor_atom_obj: # Should not happen if all_atoms_dict is complete
							continue 
							
						# 5. Precise Distance Check
						dx = target_atom_coords[0] - neighbor_atom_obj.x
						dy = target_atom_coords[1] - neighbor_atom_obj.y
						dz = target_atom_coords[2] - neighbor_atom_obj.z
						dist_sq = dx*dx + dy*dy + dz*dz # Calculate squared distance first

						if dist_sq < cutoff_radius * cutoff_radius:
							found_neighbors.append(neighbor_atom_obj)
							
	# Note: If atoms can appear in multiple cells' lists due to search span logic (less likely here),
	# converting found_neighbors (if it stored IDs) to a set then back to list of objects would ensure uniqueness.
	# Since we are adding objects and checking IDs, duplicates of the same object are not an issue unless
	# the grid itself had duplicate ID entries for some reason (which it shouldn't).
	# A final check for uniqueness based on atom_id might be warranted if strictness is needed,
	# though the current logic should generally avoid adding the same Atom object multiple times
	# unless an atom could be retrieved from different cells (which it can't, an atom belongs to one cell).
	# The search iterates cells, so an atom is processed once per cell it's in.
	# The check `atom_id == target_atom_id` is important.
	
	# To ensure unique atom *objects* if some complex scenario could lead to it:
	# unique_neighbors = {atom.id: atom for atom in found_neighbors}
	# return list(unique_neighbors.values())
	# For now, returning found_neighbors directly. If performance becomes an issue with many
	# redundant distance calculations for very dense overlapping search spheres, then optimize.
	# The current approach iterates candidate cells and then atoms in those cells.
	
	# To ensure we don't have duplicate *Atom objects* if somehow they were added multiple times
	# (e.g. if an atom was listed in multiple cells, which create_spatial_grid tries to avoid),
	# a final filter for unique objects might be useful.
	# However, each atom_id is unique in all_atoms_dict. The list `found_neighbors` contains Atom objects.
	# If the same Atom object is truly within radius and picked up from considering different
	# *neighboring_cells_of_the_target_atom's_cell*, it should only be added once due to `atom_id == target_atom_id`.
	# The critical part is that `atom_ids_in_cell` are unique per cell, and each atom is assigned to one cell.
	# The search loops through cells, then atoms in those cells. Distance check is final.
	# So, `found_neighbors` should contain unique Atom objects.

	return found_neighbors