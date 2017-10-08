#!/usr/bin/env python3
'''
Expansion of ClusterPicker (Manon Ragonnet & Emma Hodcroft):
* The algorithms implemented in ClusterPicker were at least quadratic in time
  complexity, so they have been implemented in linear time here
* ClusterPicker requires clusters to correspond to entire subtrees. We have
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
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0; descendant.num_leaves = 0
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
        child_nodes = node.child_nodes()
        assert len(child_nodes) in {0,2}, "ERROR: Multifurcating tree. Resolve polytomies first"
        if len(child_nodes) == 0:
            leaves.add(node.taxon.label)
        if node.label is None:
            node.label = 1.
        else:
            node.label = float(node.label)
            if node.label < support:
                node.edge_length = float('inf') # don't allow low-support edges
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
            node.left_dist = max(child_nodes[0].left_dist,child_nodes[0].right_dist) + child_nodes[0].edge_length
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

# min_clusters_threshold_max, but all clusters must define a subtree
def min_clusters_threshold_max_subtree(tree,threshold,support):
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
            node.left_dist = max(child_nodes[0].left_dist,child_nodes[0].right_dist) + child_nodes[0].edge_length
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

# split leaves into minimum number of clusters such that the average leaf pairwise distance is below some threshold
# I think this is incorrect because I choose to cut the child with larger total distance (not always correct)
def min_clusters_threshold_avg(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.postorder_node_iter():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        child_nodes = node.child_nodes()
        if len(child_nodes) == 0:
            node.num_leaves = 1
            node.total_dist = 0
        else:
            node.num_leaves = child_nodes[0].num_leaves + child_nodes[1].num_leaves
            if node.num_leaves == 0:
                cut(node); continue
            nl = child_nodes[0].num_leaves; nr = child_nodes[1].num_leaves
            dl = float(child_nodes[0].total_dist); dr = float(child_nodes[1].total_dist)
            el = child_nodes[0].edge_length; er = child_nodes[1].edge_length
            node.total_dist = nr*dl + nl*dr + (nl*nr)*(el+er)

            # if my kids are screwing things up, cut out the longer one
            if nl != 0 and nr != 0 and node.total_dist/(nl*nr) > threshold:
                if dl > dr:
                    node.left_dist = 0; node.num_leaves -= child_nodes[0].num_leaves
                    cluster = cut(child_nodes[0])
                else:
                    node.right_dist = 0; node.num_leaves -= child_nodes[1].num_leaves
                    cluster = cut(child_nodes[1])

                # add cluster
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# min_clusters_threshold_avg, but all clusters must define a subtree
def min_clusters_threshold_avg_subtree(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = []
    for node in tree.postorder_node_iter():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        child_nodes = node.child_nodes()
        if len(child_nodes) == 0:
            node.num_leaves = 1
            node.total_dist = 0
        else:
            node.num_leaves = child_nodes[0].num_leaves + child_nodes[1].num_leaves
            if node.num_leaves == 0:
                cut(node); continue
            nl = child_nodes[0].num_leaves; nr = child_nodes[1].num_leaves
            dl = float(child_nodes[0].total_dist); dr = float(child_nodes[1].total_dist)
            el = child_nodes[0].edge_length; er = child_nodes[1].edge_length
            node.total_dist = nr*dl + nl*dr + (nl*nr)*(el+er)

            # if my kids are screwing things up, cut out the longer one
            if nl != 0 and nr != 0 and node.total_dist/(nl*nr) > threshold:
                node.left_dist = 0; node.right_dist = 0; node.num_leaves = 0
                cluster_l = cut(child_nodes[0])
                cluster_r = cut(child_nodes[1])

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

METHODS = {'max':min_clusters_threshold_max, 'avg':min_clusters_threshold_avg, 'max_subtree':min_clusters_threshold_max_subtree, 'avg_subtree':min_clusters_threshold_avg_subtree}
if __name__ == "__main__":
    # parse user arguments
    import argparse
    import dendropy
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input Tree File")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Length Threshold")
    parser.add_argument('-s', '--support', required=False, type=float, default=0, help="Branch Support Threshold")
    parser.add_argument('-m', '--method', required=False, type=str, default='max', help="Clustering Method (options: %s)" % ', '.join(sorted(METHODS.keys())))
    args = parser.parse_args()
    if args.input == 'stdin':
        from sys import stdin; infile = stdin
    else:
        infile = open(args.input)
    trees = [dendropy.Tree.get(data=line.strip(),schema='newick',preserve_underscores=True) for line in infile.read().strip().splitlines()]
    assert args.method.lower() in METHODS, "ERROR: Invalid method: %s" % args.method
    assert args.threshold >= 0, "ERROR: Length threshold must be at least 0"
    assert args.support >= 0, "ERROR: Branch support must be at least 0"

    # run algorithm
    for t,tree in enumerate(trees):
        clusters = METHODS[args.method.lower()](tree,args.threshold,args.support)
        f = open('%s.tree%d.list.txt' % (args.input,t+1), 'w')
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