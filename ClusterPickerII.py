#!/usr/bin/env python3
from math import log
from queue import PriorityQueue,Queue

# merge two sorted lists into a sorted list
def merge_two_sorted_lists(x,y):
    out = []; i = 0; j = 0
    while i < len(x) and j < len(y):
        if x[i] < y[j]:
            out.append(x[i]); i+= 1
        else:
            out.append(y[j]); j += 1
    while i < len(x):
        out.append(x[i]); i += 1
    while j < len(y):
        out.append(y[j]); j += 1
    return out

# merge multiple sorted lists into a sorted list
def merge_multi_sorted_lists(lists):
    pq = PriorityQueue()
    for l in range(len(lists)):
        if len(lists[l]) != 0:
            pq.put((lists[l][0],l))
    inds = [1 for _ in range(len(lists))]
    out = []
    while not pq.empty():
        d,l = pq.get(); out.append(d)
        if inds[l] < len(lists[l]):
            pq.put((lists[l][inds[l]],l)); l += 1
    return out

# get the median of a sorted list
def median(x):
    if len(x) % 2 != 0:
        return x[int(len(x)/2)]
    else:
        return (x[int(len(x)/2)]+x[int(len(x)/2)-1])/2

# convert p-distance to Jukes-Cantor distance
def p_to_jc(d,seq_type):
    b = {'dna':3./4., 'protein':19./20.}[seq_type]
    return -1*b*log(1-(d/b))

# cut out the current node's subtree (by setting all nodes' DELETED to True) and return list of leaves
def cut(node):
    cluster = []
    descendants = Queue(); descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.branch_length = 0
        if descendant.is_terminal():
            cluster.append(descendant.name)
        else:
            for c in descendant.clades:
                descendants.put(c)
    return cluster

# initialize properties of input tree and return set containing taxa of leaves
def prep(tree,support):
    tree.rooted = True
    leaves = set()
    for node in tree.find_clades(order='postorder'):
        assert len(node.clades) in {0,2}, "ERROR: Multifurcating tree. Resolve polytomies first"
        node.DELETED = False
        if node.is_terminal():
            leaves.add(node.name)
        else:
            try:
                node.name = float(node.name)
            except:
                node.name = 100. # give edges without support values support 100
            if node.name < support: # don't allow low-support edges
                node.branch_length = float('inf')
    return leaves

# split leaves into minimum number of clusters such that the maximum leaf pairwise distance is below some threshold
def min_clusters_threshold_max(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.find_clades(order='postorder'):
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_terminal():
            node.left_dist = 0; node.right_dist = 0
        else:
            if node.clades[0].DELETED and node.clades[1].DELETED:
                cut(node); continue
            if node.clades[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(node.clades[0].left_dist,node.clades[0].right_dist) + node.clades[0].branch_length
            if node.clades[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(node.clades[1].left_dist,node.clades[1].right_dist) + node.clades[1].branch_length

            # if my kids are screwing things up, cut out the longer one
            if node.left_dist + node.right_dist > threshold:
                if node.left_dist > node.right_dist:
                    cluster = cut(node.clades[0])
                    node.left_dist = 0
                else:
                    cluster = cut(node.clades[1])
                    node.right_dist = 0

                # add cluster
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# median leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_med_clade(tree,threshold,support):
    leaves = prep(tree,support)
    # bottom-up traversal to compute median pairwise distances
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            node.med_pair_dist = 0
            node.leaf_dists = [0]
            node.pair_dists = []
        else:
            l_leaf_dists = [d + node.clades[0].branch_length for d in node.clades[0].leaf_dists]
            r_leaf_dists = [d + node.clades[1].branch_length for d in node.clades[1].leaf_dists]
            node.leaf_dists = merge_two_sorted_lists(l_leaf_dists,r_leaf_dists)
            node.pair_dists = merge_multi_sorted_lists([node.clades[0].pair_dists,node.clades[1].pair_dists] + [[l+r for r in r_leaf_dists] for l in l_leaf_dists])
            if node.pair_dists[-1] == float('inf'):
                node.med_pair_dist = float('inf')
            else:
                node.med_pair_dist = median(node.pair_dists)

    # top-down traversal to cut out clusters
    clusters = []
    traverse = Queue(); traverse.put(tree.root)
    while not traverse.empty():
        node = traverse.get()
        if node.med_pair_dist <= threshold:
            clusters.append(cut(node))
        else:
            traverse.put(node.clades[0]); traverse.put(node.clades[1])
    return clusters

# average leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_avg_clade(tree,threshold,support):
    leaves = prep(tree,support)
    # bottom-up traversal to compute average pairwise distances
    for node in tree.find_clades(order='postorder'):
        node.total_pair_dist = 0; node.total_leaf_dist = 0
        if node.is_terminal():
            node.num_leaves = 1
            node.avg_pair_dist = 0
        else:
            node.num_leaves = node.clades[0].num_leaves + node.clades[1].num_leaves
            node.total_pair_dist = node.clades[0].total_pair_dist + node.clades[1].total_pair_dist + (node.clades[0].total_leaf_dist*node.clades[1].num_leaves + node.clades[1].total_leaf_dist*node.clades[0].num_leaves)
            node.total_leaf_dist = (node.clades[0].total_leaf_dist + node.clades[0].branch_length*node.clades[0].num_leaves) + (node.clades[1].total_leaf_dist + node.clades[1].branch_length*node.clades[1].num_leaves)
            node.avg_pair_dist = node.total_pair_dist/((node.num_leaves*(node.num_leaves-1))/2)

    # top-down traversal to cut out clusters
    clusters = []
    traverse = Queue(); traverse.put(tree.root)
    while not traverse.empty():
        node = traverse.get()
        if node.avg_pair_dist <= threshold:
            clusters.append(cut(node))
        else:
            traverse.put(node.clades[0]); traverse.put(node.clades[1])
    return clusters

# clusters must define clades, and clades are joined if at least one pair of leaves across the two clades are within the threshold distance
def min_clusters_threshold_single_linkage_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.find_clades(order='postorder'):
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted min distance to leaf
        if node.is_terminal():
            node.left_dist = 0; node.right_dist = 0
        else:
            if node.clades[0].DELETED and node.clades[1].DELETED:
                cut(node); continue
            if node.clades[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = min(node.clades[0].left_dist,node.clades[0].right_dist) + node.clades[0].branch_length
            if node.clades[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = min(node.clades[1].left_dist,node.clades[1].right_dist) + node.clades[1].branch_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(node.clades[0])
                node.left_dist = 0
                cluster_r = cut(node.clades[1])
                node.right_dist = 0

                # add cluster
                for cluster in (cluster_l,cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.find_clades(order='postorder'):
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_terminal():
            node.left_dist = 0; node.right_dist = 0
        else:
            if node.clades[0].DELETED and node.clades[1].DELETED:
                cut(node); continue
            if node.clades[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(node.clades[0].left_dist,node.clades[0].right_dist) + node.clades[0].branch_length
            if node.clades[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(node.clades[1].left_dist,node.clades[1].right_dist) + node.clades[1].branch_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(node.clades[0])
                node.left_dist = 0
                cluster_r = cut(node.clades[1])
                node.right_dist = 0

                # add cluster
                for cluster in (cluster_l,cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

METHODS = {'max':min_clusters_threshold_max, 'max_clade':min_clusters_threshold_max_clade, 'avg_clade':min_clusters_threshold_avg_clade, 'med_clade':min_clusters_threshold_med_clade, 'single_linkage_clade':min_clusters_threshold_single_linkage_clade}
if __name__ == "__main__":
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input Tree File")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Length Threshold")
    parser.add_argument('-s', '--support', required=False, type=float, default=float('-inf'), help="Branch Support Threshold")
    parser.add_argument('-m', '--method', required=False, type=str, default='max_clade', help="Clustering Method (options: %s)" % ', '.join(sorted(METHODS.keys())))
    args = parser.parse_args()
    assert args.method.lower() in METHODS, "ERROR: Invalid method: %s" % args.method
    assert args.threshold >= 0, "ERROR: Length threshold must be at least 0"
    assert args.support >= 0 or args.support == float('-inf'), "ERROR: Branch support must be at least 0"
    if args.input == 'stdin':
        from sys import stdin; infile = stdin
    else:
        infile = open(args.input)
    from Bio import Phylo
    trees = [tree for tree in Phylo.parse(infile,'newick')]

    # run algorithm
    for t,tree in enumerate(trees):
        clusters = METHODS[args.method.lower()](tree,args.threshold,args.support)
        f = open('%s.tree_%d.method_%s.thresh_%f.list.txt' % (args.input,t+1,args.method.lower(),args.threshold), 'w')
        f.write('SequenceName\tClusterNumber\n')
        cluster_num = 1
        for cluster in clusters:
            if len(cluster) == 1:
                f.write('%s\t-1\n' % list(cluster)[0])
            else:
                for l in cluster:
                    f.write('%s\t%d\n' % (l,cluster_num))
                cluster_num += 1
        f.close()
