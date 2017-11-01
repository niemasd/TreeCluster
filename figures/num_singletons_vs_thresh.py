#! /usr/bin/env python3
'''
Niema Moshiri 2017
Generate plots of number of singletons vs. threshold for ClusterPicker-II
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
hcv = {'t':[],'n':[]}
for l in [line.strip() for line in open(d+'LANL_web_2008_HCV_core_num_singletons_vs_thresh.tsv').read().strip().splitlines()]:
    t,n = l.split('\t'); hcv['t'].append(float(t)); hcv['n'].append(int(n))
hiv_lanl = {'t':[],'n':[]}
for l in [line.strip() for line in open(d+'LANL_web_2016_HIV1_pol_ALL_num_singletons_vs_thresh.tsv').read().strip().splitlines()]:
    t,n = l.split('\t'); hiv_lanl['t'].append(float(t)); hiv_lanl['n'].append(int(n))
hiv_sd = {'t':[],'n':[]}
for l in [line.strip() for line in open(d+'SD_SubtypeB_639_num_singletons_vs_thresh.tsv').read().strip().splitlines()]:
    t,n = l.split('\t'); hiv_sd['t'].append(float(t)); hiv_sd['n'].append(int(n))
flu = {'t':[],'n':[]}
for l in [line.strip() for line in open(d+'Influenza_TypeA_Human_SegmentPB2_FullLengthOnly_num_singletons_vs_thresh.tsv').read().strip().splitlines()]:
    t,n = l.split('\t'); flu['t'].append(float(t)); flu['n'].append(int(n))

# plot number of singletons vs. threshold
fig = plt.figure()
line_width = 3
plt.plot(hcv['t'],hcv['n'],color=pal['HCV'],linewidth=line_width)
plt.plot(hiv_lanl['t'],hiv_lanl['n'],color=pal['HIV_LANL'],linewidth=line_width)
plt.plot(hiv_sd['t'],hiv_sd['n'],color=pal['HIV_SD'],linewidth=line_width)
plt.plot(flu['t'],flu['n'],color=pal['flu'],linewidth=line_width)
legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
sns.plt.xlabel('Distance Threshold (Per-Site Expected Number of Mutations)',fontsize=14)
sns.plt.ylabel('Number of Singletons',fontsize=14)
sns.plt.title('Number of Singletons vs. Distance Threshold',fontsize=18)
sns.plt.show()
fig.savefig('num_singletons_vs_thresh.pdf', format='pdf', bbox_extra_artists=(legend,), bbox_inches='tight')
plt.close()