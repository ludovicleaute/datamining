"""

Example de test de l'algo des kmeans. 

La fonction des kmeans est celle de Scipy (librairie), elle prend en input un array de coordonnées, voici un exemple:

[[  45.318  397.   ]
 [  53.543  489.   ]
 [  39.683  342.   ]
 ..., 
 [  51.644  459.   ]
 [ 103.358  915.   ]
 [  27.567  249.   ]]


Et retourne (entre autres) un array avec des données sur le clusters: array([1, 1, 1, ..., 1, 0, 1]
label[i] is the code or index of the centroid the i’th observation is closest to.


blabla

"""


from numpy import array
from scipy.cluster.vq import vq, kmeans2, whiten

file_name = "uniprot.tab"

def create_mass_length_array():

	list_tmp = []

	f = open(file_name,"r")
	headers = f.readline().rstrip().split('\t')

	for line in f.readlines():
		if not line :
			break
	    
		line_tmp = line.rstrip().split('\t')
		list_tmp.append( [ float(line_tmp[1]) ,  float(line_tmp[2]) ] )

	f.close()

	return array( list_tmp )
	

def launch_kmeans(array_mass_length):

	whitened = whiten(array_mass_length)
	book = array((whitened[0],whitened[2]))

	k = kmeans2(whitened,book)

	return k

array_mass_length = create_mass_length_array()
print launch_kmeans(array_mass_length)
