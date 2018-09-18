"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

__author__ = 'Allison MacLeay'

FILE = '/Users/Admin/Documents/mgh/projects/pipeline_files/fusions/marty/known_values.txt'

features = ['remainder_exon19_average_coverage','remainder_exon19_max_coverage', 'candidate_increased_max_coverage', 'candidate_increased_average_coverage', 'background_alk_max_coverage', 'background_alk_average_coverage']

df = pd.read_csv(FILE, sep='\t')

# Separating out the features
x = df.loc[:, features].values
# Separating out the target
y = df.loc[:,['positive']].values
# Standardizing the features
x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, df[['positive']]], axis = 1)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [1, 0]
colors = ['r', 'g'] # , 'b']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['positive'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()