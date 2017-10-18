# ClusterPicker-II
ClusterPicker-II is a tool that, given a phylogenetic tree *T* and a length threshold *L*, finds the minimum number of clusters of the leaves of *T* such that some user-specific constraint is met in each cluster. The default constraint is that the diameter of each cluster cannot exceed *L* (i.e., for each cluster, all pairs of leaves in the cluster are at most *L* apart).

ClusterPicker-II was motivated by [ClusterPicker](https://github.com/emmahodcroft/cluster-picker-and-cluster-matcher). In ClusterPicker, there is an added constraint that all clusters must define a clade, which is an available option in ClusterPicker-II as well but is not recommended based on simulation experiments (an implied constraint of this method is that trees must be rooted). Further, ClusterPicker's implementation is cubic time with respect to the number of leaves of the tree in the worst case, whereas ClusterPicker-II is linear time.

## Usage
```bash
usage: ClusterPickerII.py [-h] [-i INPUT] -t THRESHOLD [-s SUPPORT]
                          [-m METHOD]

Expansion of ClusterPicker (Manon Ragonnet & Emma Hodcroft): * The algorithms
implemented in ClusterPicker were at least quadratic in time complexity, so
they have been implemented in linear time here * ClusterPicker requires
clusters to correspond to entire clades. We have added algorithms that do not
have this restriction (also linear-time)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Tree File (default: stdin)
  -t THRESHOLD, --threshold THRESHOLD
                        Length Threshold (default: None)
  -s SUPPORT, --support SUPPORT
                        Branch Support Threshold (default: 0)
  -m METHOD, --method METHOD
                        Clustering Method (options: max, max_clade) (default:
                        max)
```

## Requirements
* [Dendropy](https://www.dendropy.org/)
