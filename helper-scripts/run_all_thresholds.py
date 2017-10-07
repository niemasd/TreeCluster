#!/usr/bin/env python3
'''
Run ClusterPicker-II on a given tree using all relevant thresholds:
* Start with threshold = 0
* Increment threshold by 0.001 (smallest branch length would take forever)
* End with threshold = longest pairwise branch length + smallest branch length
'''
INCREMENT = 0.001 # smallest bl was ~10^-9 and max distance was ~1.1, so that would take forever

# get shortest branch length and maximum pairwise distance
def check(tree):
    shortest_bl = float('inf')
    max_dist = float('-inf')
    for node in tree.postorder_node_iter():
        child_nodes = node.child_nodes()
        assert len(child_nodes) in {0,2}, "ERROR: Multifurcating tree. Resolve polytomies first"
        if node.edge_length is not None and node.edge_length > 0 and node.edge_length < shortest_bl:
            shortest_bl = node.edge_length
        if len(child_nodes) == 0:
            node.closest_leaf_dist = 0
        else:
            left_dist = child_nodes[0].closest_leaf_dist + child_nodes[0].edge_length
            right_dist = child_nodes[1].closest_leaf_dist + child_nodes[1].edge_length
            if left_dist + right_dist > max_dist:
                max_dist = left_dist + right_dist
            node.closest_leaf_dist = max(left_dist,right_dist)
    return shortest_bl,max_dist

if __name__ == "__main__":
    from os import rename
    from os.path import realpath
    from subprocess import call
    import argparse
    import dendropy
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, help="Input Tree File")
    parser.add_argument('-m', '--method', required=True, type=str, help="Clustering Method")
    args = parser.parse_args()
    tree = dendropy.Tree.get(data=open(args.input).read().strip(),schema='newick',preserve_underscores=True)
    shortest_bl,max_dist = check(tree)
    end_threshold = max_dist + shortest_bl
    threshold = 0.
    print("Running all thresholds from 0 to max pairwise distance (%f) with increment %f on tree: %s" % (end_threshold,INCREMENT,args.input))
    while threshold < end_threshold:
        print("Threshold: %f" % threshold)
        call(['%s/ClusterPickerII.py' % '/'.join(realpath(__file__).split('/')[:-2]),'-i',args.input,'-t',str(threshold),'-m',args.method])
        rename('%s.tree1.list.txt' % args.input, '%s.%f.list.txt' % (args.input,threshold))
        threshold += INCREMENT
