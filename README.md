# ClusterPicker-II
ClusterPicker-II is a tool that, given a tree *T* (Newick format) and a distance threshold *t*, finds the minimum number of clusters of the leaves of *T* such that some user-specific constraint is met in each cluster. The default constraint is that the diameter of each cluster cannot exceed *t* (i.e., for each cluster, all pairs of leaves in the cluster are at most *t* apart) and that clusters must form clades in the tree.

ClusterPicker-II was motivated by [ClusterPicker](https://github.com/emmahodcroft/cluster-picker-and-cluster-matcher).

## Usage
```bash
usage: ClusterPickerII.py [-h] [-i INPUT] -t THRESHOLD [-s SUPPORT]
                          [-m METHOD]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Tree File (default: stdin)
  -t THRESHOLD, --threshold THRESHOLD
                        Length Threshold (default: None)
  -s SUPPORT, --support SUPPORT
                        Branch Support Threshold (default: -inf)
  -m METHOD, --method METHOD
                        Clustering Method (options: max, max_clade) (default:
                        max_clade)
```

## Requirements
* [Biopython](http://biopython.org/)
