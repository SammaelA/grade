import pandas as pd 
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"
# loading dataset into Pandas DataFrame
#df = pd.read_csv(url
#                 , names=['sepal length','sepal width','petal length','petal width','target'])
#df.head()
df = pd.read_csv('test.csv')
print(df)
#features = ['sepal length', 'sepal width', 'petal length', 'petal width']
features = list(df)
#for i in range(512):
#    features.append('column_'+str(i)) 
x = df.loc[:, features].values

x = StandardScaler().fit_transform(x)
pd.DataFrame(data = x, columns = features).head()
pca = PCA(n_components=3)
principalComponents = pca.fit_transform(x)

principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2', 'principal component 3'])
finalDf = pd.concat([principalDf, df[['target']]], axis = 1)
finalDf.head(5)

fig = plt.figure(figsize = (8,8))
ax = mplot3d.Axes3D(fig)

ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)

ax.set_title('3 Component PCA', fontsize = 20)

print(finalDf.loc[:, 'principal component 1'])
targets = [0, 1, 2, 3]
colors = ['r', 'g', 'b', 'magenta']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    print(finalDf['target'])
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()

plt.show()
