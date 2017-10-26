#!/usr/bin/env python3
'''
Expansion of ClusterPicker (Manon Ragonnet & Emma Hodcroft):
* The algorithms implemented in ClusterPicker were at least quadratic in time
  complexity, so they have been implemented in linear time here
* ClusterPicker requires clusters to correspond to entire clades. We have
  added algorithms that do not have this restriction (also linear-time)
'''
from math import log
from queue import Queue

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
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0
        desc_children = descendant.child_nodes()
        if len(desc_children) == 0:
            cluster.append(descendant.taxon.label)
        else:
            for c in desc_children:
                descendants.put(c)
    return cluster

# initialize properties of input tree and return set containing taxa of leaves
def prep(tree,support):
    tree.seed_node.edge_length = 0
    leaves = set()
    for node in tree.postorder_node_iter():
        node.DELETED = False
        if node.is_leaf():
            leaves.add(node.taxon.label)
        if node.label is None:
            node.label = 100.
        else:
            try:
                node.label = float(node.label)
                if node.label < support:
                    node.edge_length = float('inf') # don't allow low-support edges
            except:
                node.label = 100. # ignore internal node labels if they're not support values
    return leaves

# split leaves into minimum number of clusters such that the maximum leaf pairwise distance is below some threshold
def min_clusters_threshold_max(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.postorder_node_iter():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        child_nodes = node.child_nodes()
        if len(child_nodes) == 0:
            node.left_dist = 0; node.right_dist = 0
        else:
            if child_nodes[0].DELETED and child_nodes[1].DELETED:
                cut(node)
            if child_nodes[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(child_nodes[0].left_dist,child_nodes[0].right_dist) + child_nodes[0].edge_length
            if child_nodes[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(child_nodes[1].left_dist,child_nodes[1].right_dist) + child_nodes[1].edge_length

            # if my kids are screwing things up, cut out the longer one
            if node.left_dist + node.right_dist > threshold:
                if node.left_dist > node.right_dist:
                    cluster = cut(child_nodes[0])
                    node.left_dist = 0
                else:
                    cluster = cut(child_nodes[1])
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

# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.postorder_node_iter():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        child_nodes = node.child_nodes()
        if len(child_nodes) == 0:
            node.left_dist = 0; node.right_dist = 0
        else:
            if child_nodes[0].DELETED and child_nodes[1].DELETED:
                cut(node)
            if child_nodes[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(child_nodes[0].left_dist,child_nodes[0].right_dist) + child_nodes[0].edge_length
            if child_nodes[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(child_nodes[1].left_dist,child_nodes[1].right_dist) + child_nodes[1].edge_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(child_nodes[0])
                node.left_dist = 0
                cluster_r = cut(child_nodes[1])
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

METHODS = {'max':min_clusters_threshold_max, 'max_clade':min_clusters_threshold_max_clade}
if __name__ == "__main__":
    # parse user arguments
    import argparse
    import dendropy
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input Tree File")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Length Threshold")
    parser.add_argument('-s', '--support', required=False, type=float, default=0, help="Branch Support Threshold")
    parser.add_argument('-m', '--method', required=False, type=str, default='max_clade', help="Clustering Method (options: %s)" % ', '.join(sorted(METHODS.keys())))
    args = parser.parse_args()
    assert args.method.lower() in METHODS, "ERROR: Invalid method: %s" % args.method
    assert args.threshold >= 0, "ERROR: Length threshold must be at least 0"
    assert args.support >= 0, "ERROR: Branch support must be at least 0"
    if args.input == 'stdin':
        from sys import stdin; infile = stdin
    else:
        infile = open(args.input)
    trees = [dendropy.Tree.get(data=line.strip(),schema='newick',preserve_underscores=True) for line in infile.read().strip().splitlines()]
    for tree in trees:
        tree.suppress_unifurcations()
        tree.resolve_polytomies()

    # run algorithm
    for t,tree in enumerate(trees):
        clusters = METHODS[args.method.lower()](tree,args.threshold,args.support)
        f = open('%s.tree%d.thresh%f.list.txt' % (args.input,t+1,args.threshold), 'w')
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
