from scipy.spatial import ConvexHull
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import numpy as np
from itertools import product
from collections import defaultdict
import multiprocessing as mp
from sys import argv

def worker(info):
    i, actives = info
    if i == 0 or i == (1 << 16) - 1:
        return []
    active_vertices = np.array(actives).reshape(2,2,2,2)
    corner_vertices = np.array(np.nonzero(active_vertices)).T
    midpoint_vertices = set()
    for vertex in corner_vertices:
        for j in range(4):
            copy_vert = vertex.copy()
            copy_vert[j] = 1 - copy_vert[j]
            if not active_vertices[tuple(copy_vert)]:
                midpoint_vertices.add(tuple((copy_vert + vertex) / 2))
    midpoint_vertices = np.array(list(midpoint_vertices))
    all_vertices = np.vstack((corner_vertices, midpoint_vertices))
    hull = ConvexHull(all_vertices)
    hyperplanes = []
    idx_map = {}
    hyperplane_map = {}
    for j, equation in enumerate(hull.equations):
        split_eq = (equation[:-1], equation[-1])
        if not (np.count_nonzero(split_eq[0]) == 1 and split_eq[1] in [0, -1]):
            if tuple(equation) in hyperplane_map:
                idx_map[j] = hyperplane_map[tuple(equation)]
            else:
                hyperplane_map[tuple(equation)] = idx_map[j] = len(hyperplanes)
                hyperplanes.append(split_eq)
    connectivity = csr_matrix((len(hyperplanes), len(hyperplanes)), dtype=np.int8)
    for j, neighbors in enumerate(hull.neighbors):
        if j in idx_map:
            p = idx_map[j]
            for neighbor in neighbors:
                if neighbor in idx_map:
                    q = idx_map[neighbor]
                    connectivity[p, q] = 1
                    connectivity[q, p] = 1
    _, labels = connected_components(connectivity)
    labelled = defaultdict(list)
    for label, hyperplane in zip(labels, hyperplanes):
        labelled[label].append(hyperplane)
    return list(labelled.values())

def write_table(table):
    with open(argv[1], 'w') as f:
        f.write('#include "stmesh/table.hpp"\n\n')
        f.write('namespace stmesh::detail {\n')
        f.write('const double kHyperplanes[] = {\n')
        hyperplane_idx = 0
        hyperplane_idxs = []
        component_idx = 0
        component_idxs = []
        for components in table:
            component_idxs.append(component_idx)
            component_idx += len(components)
            for i, component in enumerate(components):
                hyperplane_idxs.append(hyperplane_idx)
                hyperplane_idx += len(component) * 5
                for hyperplane in component:
                    f.write(f'  {hyperplane[0][0]}, {hyperplane[0][1]}, {hyperplane[0][2]}, {hyperplane[0][3]}, {hyperplane[1]}, \n')
                if i == len(components) - 1:
                    f.write('  /*################################################################*/\n')
                else:
                    f.write('  /******************************************************************/\n')
        component_idxs.append(component_idx)
        hyperplane_idxs.append(hyperplane_idx)
        f.write('};\n\n')
        f.write('const unsigned kHyperplaneIdxs[] = {\n')
        for i in hyperplane_idxs:
            f.write(f'  {i},\n')
        f.write('};\n\n')
        f.write('const unsigned kComponentIdxs[] = {\n')
        for i in component_idxs:
            f.write(f'  {i},\n')
        f.write('};\n')
        f.write('} // namespace stmesh::detail\n')

if __name__ == '__main__':
    with mp.Pool(mp.cpu_count()) as pool:
        table = pool.map(worker, enumerate(product(*(([False, True],)*16))))
    write_table(table)
