#! /usr/bin/env python3
'''
Niema Moshiri 2017
Generate plots of runtime vs. n for ClusterPicker and ClusterPicker-II
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
pal = {'cp':'#FF0000','cp2':'#0000FF'}
handles = [Patch(color=pal['cp'],label='ClusterPicker'), Patch(color=pal['cp2'],label='ClusterPicker-II')]

# data
time_cp = {
    'n': [100]*10 + [200]*10 + [500]*10 + [1000]*10 + [2000]*10 + [5000]*10,
    't': [0.765,0.913,0.833,0.841,0.811,1.043,0.872,0.72,0.898,1.27] +                              # n = 100
         [4.572,2.265,2.368,3.11,2.494,4.272,2.828,2.12,2.914,2.58] +                               # n = 200
         [20.42,24.538,23.503,22.623,19.922,10.507,14.237,12.197,6.297,11.544] +                    # n = 500
         [80.233,65.727,91.069,85.733,125.850,63.788,83.152,77.214,52.553,73.051] +                 # n = 1000
         [497.015,678.894,603.918,622.179,437.459,812.144,430.469,402.385,395.973,663.566] +        # n = 2000
         [2856.56,3232.143,3775.304,5165.85,4150.273,4683.825,4212.49,3420.913,2986.446,5754.584] + # n = 5000
         []
}
time_cp2 = {
    'n': [100]*10 + [200]*10 + [500]*10 + [1000]*10 + [2000]*10 + [5000]*10,
    't': [0.281,0.28,0.276,0.277,0.276,0.276,0.276,0.276,0.276,0.276] +  # n = 100
         [0.296,0.295,0.294,0.295,0.298,0.295,0.295,0.293,0.295,0.295] + # n = 200
         [0.349,0.35,0.35,0.352,0.35,0.349,0.349,0.349,0.351,0.351] +    # n = 500
         [0.442,0.439,0.441,0.444,0.444,0.456,0.446,0.444,0.443,0.439] + # n = 1000
         [0.625,0.629,0.62,0.62,0.619,0.625,0.622,0.624,0.625,0.625] +   # n = 2000
         [1.164,1.171,1.173,1.168,1.165,1.17,1.174,1.176,1.183,1.185] +  # n = 5000
         []
}

# plot time vs. n
fig = plt.figure()
ax = sns.pointplot(x='n',y='t',data=pd.DataFrame(time_cp),color=pal['cp'],markers=[''])
sns.pointplot(x='n',y='t',data=pd.DataFrame(time_cp2),color=pal['cp2'],markers=[''])
ax.set_yscale('log')
legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
sns.plt.xlabel('Number of Taxa',fontsize=14)
sns.plt.ylabel('Execution Time (s)',fontsize=14)
sns.plt.title('Execution Time (s) vs. Number of Taxa',fontsize=18)
sns.plt.show()
fig.savefig('time_vs_n.pdf', format='pdf', bbox_extra_artists=(legend,), bbox_inches='tight')
plt.close()