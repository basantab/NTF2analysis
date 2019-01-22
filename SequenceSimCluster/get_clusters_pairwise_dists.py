import numpy as np
from scipy import cluster
import sys

fname_dist = sys.argv[1]
fname_dec = sys.argv[2]
maxclust_n = int(sys.argv[3])


array = np.genfromtxt(fname_dist,skip_header=1)
square_Array = np.delete(array,0,1)

#data = np.fromfile(fname_dist,sep=' ')

#size = int(np.sqrt(len(data)))

#square_Array = np.reshape(data,(size,size))

triu = np.triu_indices_from(square_Array,k=1)

condensed_upper = square_Array[triu]

condensed_upper_format = [i for i in condensed_upper]

linkage_mtx = cluster.hierarchy.linkage(condensed_upper_format,method='average')

flat_clusters = cluster.hierarchy.fcluster(linkage_mtx,maxclust_n,criterion='maxclust')

flat_clusters.tofile('CLUSTERS.raw',sep=' ')

decoys = [i[:-1] for i in open(fname_dec,'r').readlines()]

model_to_cluster_file = open('model_to_cluster.list','w')
for n,i in enumerate(flat_clusters):
	model_to_cluster_file.write('%d %s\n'%(i,decoys[n]))

model_to_cluster_file.close()
