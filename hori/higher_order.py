from multiprocessing import Pool
from tqdm import tqdm

def compute_higher_order_group(args):
	"""
	Worker function to compute higher-order groups for a given group and order.

	In this adjacency-based version, we assume:
	  args = (group, adjacency, all_residues)

	'group' is a known (k-)clique, adjacency is a dict
	{res: set_of_neighbors}, and all_residues is unused except for
	consistency with the original signature.
	"""
	group, adjacency, all_residues = args

	# Intersection of neighbors of all residues in 'group'
	common_neighbors = adjacency[group[0]].copy()
	for residue in group[1:]:
		common_neighbors &= adjacency[residue]
	# Exclude members already in the group
	common_neighbors -= set(group)

	# Each candidate is adjacent to all in 'group', so form new (k+1)-cliques
	new_groups = []
	for candidate in common_neighbors:
		new_group = tuple(sorted(group + (candidate,)))
		new_groups.append(new_group)

	return new_groups


def parallel_find_higher_order_interactions(residues, residue_interactions, max_order=None, num_cores=8, parallel=False):
	"""
	Parallelized (or single-threaded) function to find higher-order residue interactions,
	with duplicate prevention.

	If max_order='converge', the function will run until no groups or only one group is found.
	If parallel=False, it runs single-threaded. If parallel=True, it uses multiprocessing.

	Returns a dictionary of the form:
	   {2: set_of_2_cliques, 3: set_of_3_cliques, 4: set_of_4_cliques, ... }
	   or, if max_order='converge', it returns up to the order where expansion converges.
	"""

	# 1) Decide how far we expand
	max_order = max_order if max_order is not None else 'converge'

	# 2) Initialize the container
	if max_order != 'converge':
		# Pre-fill up to max_order
		highest_num = int(max_order)
		higher_order_interactions = {n: set() for n in range(2, highest_num + 1)}
	else:
		higher_order_interactions = {}

	# 3) Build adjacency and collect all 2-cliques (edges)
	adjacency = {}
	edge_set = set()  # This will store 2-cliques as (min_node, max_node)
	for (r1, r2) in residue_interactions.keys():
		adjacency.setdefault(r1, set()).add(r2)
		adjacency.setdefault(r2, set()).add(r1)
		edge_set.add(tuple(sorted([r1, r2])))

	all_residues = set(residues.keys())

	# Place 2-cliques in the dictionary so we can expand them
	if 2 in higher_order_interactions:
		higher_order_interactions[2] = edge_set
	else:
		# If in 'converge' mode, add them dynamically
		if max_order == 'converge':
			higher_order_interactions[2] = edge_set

	# 4) Expand cliques from k=2 upward
	order = 3
	while max_order == 'converge' or order <= int(max_order):
		# If we have a known (order-1)-clique set, expand from it
		if (order - 1) in higher_order_interactions and higher_order_interactions[order - 1]:
			prev_order_groups = higher_order_interactions[order - 1]
		else:
			# Fallback to 2-cliques if needed
			prev_order_groups = edge_set

		print(f"Processing {order}-residue interactions ({len(prev_order_groups)} groups)...")

		# Prepare tasks for either parallel or single-thread
		task_args = [(group, adjacency, all_residues) for group in prev_order_groups]

		if parallel:
			# ----------- MULTIPROCESSING -----------
			with Pool(processes=num_cores) as pool:
				results = list(
					tqdm(
						pool.imap_unordered(compute_higher_order_group, task_args),
						total=len(task_args),
						desc=f"{order}-Residue Progress"
					)
				)
		else:
			# ----------- SINGLE-THREAD -----------
			results = []
			for args in tqdm(task_args, desc=f"{order}-Residue Progress", total=len(task_args)):
				results.append(compute_higher_order_group(args))

		# Flatten the result lists and deduplicate
		current_order_groups = set()
		for sublist in results:
			current_order_groups.update(sublist)

		# Store in dictionary
		if max_order == 'converge':
			higher_order_interactions[order] = current_order_groups
			# Convergence check
			if len(current_order_groups) <= 1:
				print(
					f"Convergence achieved at {order}-residue interactions "
					f"with {len(current_order_groups)} group(s)."
				)
				break
		else:
			higher_order_interactions[order] = current_order_groups

		print(f"Completed {order}-residue interactions: {len(current_order_groups)} unique groups\n")

		# Stop if we've hit a numeric max
		if max_order != 'converge' and order == int(max_order):
			break

		order += 1

	return higher_order_interactions
