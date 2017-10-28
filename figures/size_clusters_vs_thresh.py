#! /usr/bin/env python3
'''
Niema Moshiri 2017
Generate plots of number of clusters vs. threshold for ClusterPicker-II
'''
# imports
from matplotlib import rcParams
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# settings
sns.set_style("ticks")
rcParams['font.family'] = 'serif'
pal = {'HCV':'#FF0000','HIV_LANL':'#00FF00','HIV_SD':'#0000FF','flu':'#FFFF00'}
handles = [Patch(color=pal['HCV'],label='Hepatitis C (HCV)'), Patch(color=pal['HIV_LANL'],label='HIV (Los Alamos)'), Patch(color=pal['HIV_SD'],label='HIV (San Diego)'), Patch(color=pal['flu'],label='Influenza')]

# load data
d = '/'.join(__file__.split('/')[:-2]) + '/real_datasets/'
hcv = {'t':[],'s':[]}
for l in [line.strip() for line in open(d+'LANL_web_2008_HCV_core_size_clusters_vs_thresh.tsv').read().strip().splitlines()]:
    t,s = l.split('\t'); s = [float(i) for i in s.split(',')]; hcv['t'] += [float(t)]; hcv['s'] += [sum(s)/float(len(s))]
hiv_lanl = {'t':[],'s':[]}
for l in [line.strip() for line in open(d+'LANL_web_2016_HIV1_pol_ALL_size_clusters_vs_thresh.tsv').read().strip().splitlines()]:
    t,s = l.split('\t'); s = [float(i) for i in s.split(',')]; hiv_lanl['t'] += [float(t)]; hiv_lanl['s'] += [sum(s)/float(len(s))]
hiv_sd = {'t':[],'s':[]}
for l in [line.strip() for line in open(d+'SD_SubtypeB_639_size_clusters_vs_thresh.tsv').read().strip().splitlines()]:
    t,s = l.split('\t'); s = [float(i) for i in s.split(',')]; hiv_sd['t'] += [float(t)]; hiv_sd['s'] += [sum(s)/float(len(s))]
flu = {'t':[],'s':[]}
for l in [line.strip() for line in open(d+'Influenza_TypeA_Human_SegmentPB2_FullLengthOnly_size_clusters_vs_thresh.tsv').read().strip().splitlines()]:
    t,s = l.split('\t'); s = [float(i) for i in s.split(',')]; flu['t'] += [float(t)]; flu['s'] += [sum(s)/float(len(s))]

# plot number of clusters vs. threshold
fig = plt.figure()
line_width = 3
plt.plot(hcv['t'],hcv['s'],color=pal['HCV'],linewidth=line_width)
plt.plot(hiv_lanl['t'],hiv_lanl['s'],color=pal['HIV_LANL'],linewidth=line_width)
plt.plot(hiv_sd['t'],hiv_sd['s'],color=pal['HIV_SD'],linewidth=line_width)
plt.plot(flu['t'],flu['s'],color=pal['flu'],linewidth=line_width)
legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
sns.plt.xlabel('Distance Threshold (Per-Site Expected Number of Mutations)',fontsize=14)
sns.plt.ylabel('Average Cluster Size',fontsize=14)
sns.plt.title('Average Cluster Size vs. Distance Threshold',fontsize=18)
sns.plt.show()
fig.savefig('size_clusters_vs_thresh.pdf', format='pdf', bbox_extra_artists=(legend,), bbox_inches='tight')
plt.close()