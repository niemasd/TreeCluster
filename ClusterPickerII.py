#!/usr/bin/env python3
'''
Expansion of ClusterPicker (Manon Ragonnet & Emma Hodcroft):
* The algorithms implemented in ClusterPicker were at least quadratic in time
  complexity, so they have been implemented in linear time here
* ClusterPicker requires clusters to correspond to entire subtrees. We have
  added algorithms that do not have this restriction (also linear-time)
'''
from queue import Queue

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
def prep(tree):
    tree.seed_node.edge_length = 0
    leaves = set()
    for node in tree.postorder_node_iter():
        node.DELETED = False
        child_nodes = node.child_nodes()
        assert len(child_nodes) in {0,2}, "ERROR: Multifurcating tree. Resolve polytomies first"
        if len(child_nodes) == 0:
            leaves.add(node.taxon.label)
    return leaves

# split leaves into minimum number of clusters such that the maximum leaf pairwise distance is below some threshold
def min_clusters_threshold_max(tree,threshold):
    leaves = prep(tree)
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

# split leaves into minimum number of clusters such that the average leaf pairwise distance is below some threshold
def min_clusters_threshold_avg(tree,threshold):
    leaves = prep(tree)
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
            nl = child_nodes[0].num_leaves; nr = child_nodes[1].num_leaves
            dl = float(child_nodes[0].total_dist); dr = float(child_nodes[1].total_dist)
            el = child_nodes[0].edge_length; er = child_nodes[1].edge_length
            node.total_dist = nr*dl + nl*dr + (nl*nr)*(el+er)

            # if my kids are screwing things up, cut out the longer one
            if node.total_dist/(nl*nr) > threshold:
                if (dl/nl) > (dr/nr):
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

METHODS = {'max':min_clusters_threshold_max, 'avg':min_clusters_threshold_avg}
if __name__ == "__main__":
    # parse user arguments
    from sys import stdin
    import argparse
    import dendropy
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--file', required=False, type=argparse.FileType('r'), default=stdin, help="Input Tree File (default: stdin)")
    parser.add_argument('-t', '--threshold', required=True, type=float, help="Length Threshold")
    parser.add_argument('-m', '--method', required=False, type=str, default='max', help="Clustering Method (options: %s)" % ', '.join(sorted(METHODS.keys())))
    args = parser.parse_args()
    trees = [dendropy.Tree.get(data=line.strip(),schema='newick') for line in args.file.read().strip().splitlines()]
    assert args.method.lower() in METHODS, "ERROR: Invalid method: %s" % args.method

    # run algorithm
    for tree in trees:
        clusters = METHODS[args.method.lower()](tree,args.threshold)
        for cluster in clusters:
            print(cluster)
        print()