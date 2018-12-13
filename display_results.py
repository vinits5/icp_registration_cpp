import os
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

matrix = np.loadtxt('transformation.txt',delimiter=',')
S_given = np.loadtxt(os.path.join(os.getcwd(),'armadillo20deg','S.txt'), delimiter = ',')
V = np.loadtxt(os.path.join(os.getcwd(),'armadillo20deg','Msampled.txt'), delimiter = ',')

S = S_given[0:500,:]

rot = matrix[0:3,0:3]
tra = np.reshape(np.transpose(matrix[0:3,3]), (3,1))
transformed = np.matmul(rot, np.transpose(S)) + tra
fig = plt.figure()
# ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')
ax.set_label("x - axis")
ax.set_label("y - axis")
ax.set_label("z - axis")

ax.plot(V[:,0], V[:,1], V[:,2], "o", color="red", ms=6, mew=0.5, label='Model Data')								# Model Data.
ax.plot(S[:,0], S[:,1], S[:,2], "o", color="blue", ms=6, mew=0, label='Sensor Data')   								# Sensor Data.
ax.plot(transformed[0,:], transformed[1,:], transformed[2,:], "o", color="green", ms=6, mew=0, label='Transformed Sensor Data')	# Transformed Sensor Data.
ax.legend(fontsize=17,markerscale=3)
plt.show()