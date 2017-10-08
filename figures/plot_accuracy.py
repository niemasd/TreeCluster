#! /usr/bin/env python3
from matplotlib import rcParams
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# settings
sns.set_style("ticks")
rcParams['font.family'] = 'serif'
pal = {'max':'#0000FF', 'max_subtree':'#FF0000'}
handles = [Patch(color=pal['max'],label='Max'), Patch(color=pal['max_subtree'],label='Max Subtree')]

# moving average
def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

# read data (CSV files should have columns: Rep, TP, TN, FP, FN, Precision, Recall)
data_max = {}
for line in open('accuracy.8.clusters.max.csv').read().strip().splitlines()[1:]:
    parts = line.split(',')
    if parts[0] not in data_max:
        data_max[parts[0]] = {'Precision':[],'Recall':[]}
    data_max[parts[0]]['Precision'].append(float(parts[-2]))
    data_max[parts[0]]['Recall'].append(float(parts[-1]))
for rep in data_max:
    data_max[rep]['Recall'],data_max[rep]['Precision'] = zip(*sorted(zip(data_max[rep]['Recall'], data_max[rep]['Precision'])))

# plot
fig = plt.figure()
axes = plt.gca()
axes.set_xlim([0,1])
axes.set_ylim([0,1])
x = []; y=[]
for rep in data_max:
    x += data_max[rep]['Recall']
    y += data_max[rep]['Precision']
x,y = zip(*sorted(zip(x,y)))
y_avg = movingaverage(y, 50)
plt.plot(x,y_avg,color=pal['max'])
legend = plt.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=True)
sns.plt.xlabel('Recall')
sns.plt.ylabel('Precision')
sns.plt.title('Precision vs. Recall (8 Clusters)')
plt.show()
