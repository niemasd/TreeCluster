#!/usr/bin/env python3
'''
Score a given query clustering against a given true clustering, where both clusterings are in the TreeCluster format.

    * AMI = Adjusted Mutual Information
    * ARI = Adjusted Rand Index
    * COM = Completeness Score
    * FMI = Fowlkes-Mallows Index
    * HCV = Compute Homogenity, Completeness, and V-Measure together
    * HOM = Homogeneity Score
    * MI  = Mutual Information
    * NMI = Normalized Mutual Information
    * VM  = V-Measure

See scikit-learn documentation for details: https://scikit-learn.org/stable/modules/classes.html#clustering-metrics
'''
from sys import stderr
try:
    from sklearn.metrics.cluster import adjusted_mutual_info_score,adjusted_rand_score,completeness_score,fowlkes_mallows_score,homogeneity_completeness_v_measure,homogeneity_score,mutual_info_score,normalized_mutual_info_score,v_measure_score
except:
    assert False, "ERROR: Unable to import sklearn. Install with: pip install scikit-learn"
METRICS = {'AMI':adjusted_mutual_info_score, 'ARI':adjusted_rand_score, 'COM':completeness_score, 'FMI':fowlkes_mallows_score, 'HCV':homogeneity_completeness_v_measure, 'HOM':homogeneity_score, 'MI':mutual_info_score, 'NMI':normalized_mutual_info_score, 'VM':v_measure_score}

# load a Cluster Picker format clustering file
def load_clusters(f):
    node_to_cluster = {}
    for line in f:
        if 'SequenceName' in line:
            continue
        n,c = [e.strip() for e in line.split()]
        node_to_cluster[n] = int(c)
    c = max(max(node_to_cluster.values())+1,1)
    for n in node_to_cluster:
        if node_to_cluster[n] == -1:
            node_to_cluster[n] = c; c += 1
    return node_to_cluster,c

if __name__ == "__main__":
    # parse args
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-q', '--query', required=True, type=argparse.FileType('r'), help="Query Clustering File")
    parser.add_argument('-r', '--reference', required=True, type=argparse.FileType('r'), help="Reference Clustering File")
    parser.add_argument('-m', '--metric', required=True, type=str, help="Scoring Method (options: %s)" % ', '.join(sorted(METRICS.keys())))
    parser.add_argument('-ns', '--no_singletons', action='store_true', help="Exclude True Singletons from Calculation")
    args,unknown = parser.parse_known_args()
    args.metric = args.metric.strip().upper()
    assert args.metric in METRICS, "ERROR: Invalid metric: %s (options: %s)" % (args.metric, ', '.join(sorted(METRICS.keys())))

    # load clusterings
    q_node_to_cluster,q_c = load_clusters(args.query)
    r_node_to_cluster,r_c = load_clusters(args.reference)
    nodes = list(r_node_to_cluster.keys())
    if args.no_singletons: # remove true singletons from nodes
        cluster_size = {}
        for n in r_node_to_cluster:
            c = r_node_to_cluster[n]
            if c not in cluster_size:
                cluster_size[c] = 0
            cluster_size[c] += 1
        nodes = [n for n in nodes if cluster_size[r_node_to_cluster[n]] > 1]
    r_cluster_list = [r_node_to_cluster[n] for n in nodes]
    q_cluster_list = []
    for n in nodes:
        if n in q_node_to_cluster:
            q_cluster_list.append(q_node_to_cluster[n])
        else: # tn93 doesn't output singletons
            q_cluster_list.append(q_c); q_c += 1

    # compute and output score
    if args.metric == 'HCV':
        h,c,v = METRICS[args.metric](r_cluster_list,q_cluster_list)
        print("HOM: %f" % h); print("COM: %f" % c); print("VM: %f" % v)
    else:
        print("%s: %f" % (args.metric, METRICS[args.metric](r_cluster_list,q_cluster_list)))
