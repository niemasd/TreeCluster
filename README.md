# TreeCluster
TreeCluster is a tool that, given a tree *T* (Newick format) and a distance threshold *t*, finds the minimum number of clusters of the leaves of *T* such that some user-specified constraint is met in each cluster. The user can also specify a branch support threshold *s* such that no pair of leaves in any cluster can be connected by branches with support below *s*.

TreeCluster was motivated by [Cluster Picker](https://github.com/emmahodcroft/cluster-picker-and-cluster-matcher).

The default method is "Max Clade" (see [Clustering Methods](#clustering-methods)). There is no explicit default distance threshold, but because Cluster Picker recommends a distance threshold of 0.045 and because the same objective function is optimized by both Cluster Picker and TreeCluster "Max Clade", we currently recommend 0.045 as well.

Note that TreeCluster can run within seconds even on ultra-large datasets, so it may make sense to use a range of thresholds and determine the appropriate choice based on the results. We intend to develop non-parametric methods of TreeCluster clustering in the future.

## Usage
```bash
usage: TreeCluster.py [-h] [-i INPUT] -t THRESHOLD [-s SUPPORT] [-m METHOD] [-tf THRESHOLD_FREE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Tree File (default: stdin)
  -t THRESHOLD, --threshold THRESHOLD
                        Length Threshold (default: None)
  -s SUPPORT, --support SUPPORT
                        Branch Support Threshold (default: -inf)
  -m METHOD, --method METHOD
                        Clustering Method (options: avg_clade, length,
                        length_clade, max, max_clade, med_clade, root_dist,
                        single_linkage_clade) (default: max_clade)
  -tf THRESHOLD_FREE, --threshold_free THRESHOLD_FREE
                        Threshold-Free Approach (options: argmax_clusters)
                        (default: None)
```

## Clustering Methods
* **Avg Clade:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The average pairwise distance between leaves in the cluster is below *t*
    2. Leaves cannot be connected by branches with support below *s*
    3. The leaves in the cluster must define a clade in *T*
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Length:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The cluster does not contain any edges above length *t*
    2. Leaves cannot be connected by branches with support below *s*
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Length Clade:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The cluster does not contain any edges above length *t*
    2. Leaves cannot be connected by branches with support below *s*
    3. The leaves in the cluster must define a clade in *T*
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Max:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The maximum pairwise distance between leaves in the cluster is below *t*
    2. Leaves cannot be connected by branches with support below *s*
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Max Clade:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The maximum pairwise distance between leaves in the cluster is below *t*
    2. Leaves cannot be connected by branches with support below *s*
    3. The leaves in the cluster must define a clade in *T*
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Med Clade:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The median pairwise distance between leaves in the cluster is below *t*
    2. Leaves cannot be connected by branches with support below *s*
    3. The leaves in the cluster must define a clade in *T*
    * For a tree with *n* leaves, this algorithm is O(*n*² log *n*) in the worst cas

* **Root Dist:** Cluster the leaves by cutting the tree at *t* distance away from the root
    * Branches with support below *s* are simply treated as infinitely long
    * For a tree with *n* leaves, this algorithm is O(*n*)

* **Single Linkage Clade:** Cluster the leaves such that the following conditions hold for each cluster:
    1. The leaves in the cluster must define a clade in *T*
    2. For all internal nodes *u* in the clade defined by the cluster, a leaf in the left subclade of *u* must be within *t* distance of a leaf in the right subclade of *u*
    3. Leaves cannot be connected by branches with support below *s*
    * For a tree with *n* leaves, this algorithm is O(*n*)

## Threshold-Free Approaches
* **Argmax Clusters:** Choose the threshold that maximizes the number of non-singleton clusters over all thresholds from 0 to *t*
    * Currently, for the sake of speed, only every 0.0001 threshold is tested (i.e., 0, 0.001, 0.002, ..., *t*)

## Requirements
* [TreeSwift](https://github.com/niemasd/TreeSwift)
