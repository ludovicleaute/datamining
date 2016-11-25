#coding:utf-8

from numpy import array
from scipy.cluster.vq import vq, kmeans2, whiten

"""

- faire structure de données avec tout (done)
- faire tableau 2D avec nom prot et [taille+longueur] (done)
- calculer moyenne
- faire une liste de cluster avec prot à l'intérieur 
- dictionnaire: pour un cluster, liste de protéine? 


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


data = {}
go_list = [] #liste de tous les goterms existants dans le jeu de donnees
k = 300 # number of cluster of kmeans method 

file_name = "uniprot.tab"


# data will associate to each protein, its properties (Mass, Length, GO terms...)
def create_data():
	
	fic = open(file_name,"r")

	headers = fic.readline().rstrip().split('\t') #noms des colonnes

	for line in fic.readlines():
		if not line:
			break

		prot = line.rstrip().split('\t')
		data[prot[0]] = {}

		for nb_col in range(1,len(prot)):

            # Mass
			if headers[nb_col] == "Mass":
				gos = prot[nb_col]
				gos = gos.replace(",",".")
				data[prot[0]]["Mass"] = float(gos)


            # Length
			if headers[nb_col] == "Length":
				gos = prot[nb_col]
				gos = gos.replace(",",".")
				data[prot[0]]["Length"] = float(gos)

            # Gene ontology (cellular component)
			if headers[nb_col] == "Gene ontology (cellular component)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GoCellularComponent"] = gos
				for go in gos:
					if go not in go_list:
						go_list.append(go)


            # Gene ontology (molecular function)
			if headers[nb_col] == "Gene ontology (molecular function)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GOMolecularFunction"] = gos
				for go in gos:
					if go not in go_list:
						go_list.append(go)
     

            # Gene ontology (biological process)
			if headers[nb_col] == "Gene ontology (biological process)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GoBiologicalProcess"] = gos
				for go in gos:
					if go not in go_list:
						go_list.append(go)


# Mass_length_array is an array containing coordinates of each prot (as x = Mass and y = Length)
def create_mass_length_array():

	list_tmp = []

	for prot in data:	    
		x = data[prot]['Length']
		y = data[prot]['Mass']
		list_tmp.append( [ x ,  y ] )

	return array( list_tmp )
	

def launch_kmeans(array_mass_length):

	# Whithened normalize data...je sais pas ce que ça veut dire mais il faut le faire
	whitened = whiten(array_mass_length)
	kmean = kmeans2(whitened,k)

	return kmean

# This table will allow to retrieve the protein associated to each point of the kmeans method
def create_2D_table():
	table = []
	
	i = 0

	for prot in data:

		table.append([])

		mass = data[prot]['Mass']
		length = data[prot]['Length']

		table[i].append(prot)
		table[i].append( [ mass , length ] )

		i += 1

	return table


def retrieve_protein(table, clusters_with_nb):
	clusters_with_name = [] 

	for i in range(k):
		clusters_with_name.append([])

	for i in range(len(clusters_with_nb)):
		protein_name = table[i][0]
		clusters_with_name[ clusters_with_nb[i] ].append(protein_name)

	return clusters_with_name






create_data()

"""
data = 

{ 
  
  PROT1 : { 'Length': 397,
            'Mass': 45.318,
            'GoBiologicalProcess': ['fructose 2,6-bisphosphate metabolic process [GO:0006003]', 'fructose metabolic process [GO:0006000]', 'regulation of glycolytic process [GO:0006110]'], 
            'GoCellularComponent': ['cytoplasm [GO:0005737]'], 
            'GOMolecularFunction': ['6-phosphofructo-2-kinase activity [GO:0003873]', 'ATP binding [GO:0005524]']   ,

  PROT2 : { ... }
}
"""

array_mass_length = create_mass_length_array()

"""
array_mass_length = 

[	 [ 724.      83.48 ]
	 [ 122.      13.994]
	 [ 112.      13.062]
	 ..., 
	 [ 132.      14.669]
	 [ 589.      68.699]
	 [ 589.      68.699]	]

"""

table = create_2D_table()


"""
table = 

[
 
 	['P53538', [23.469, 206.0]], 
 	['P26188', [21.499, 188.0]], 
 	['Q6B0W0', [12.948, 111.0]],
 	...

]
"""

kmean = launch_kmeans(array_mass_length)

centroids = kmean[0]

"""
centroids = 

[[ 0.53208188  0.53102884]
 [ 3.6624955   3.65782323]
 [ 1.57007574  1.56583634]]

"""

clusters_with_nb = kmean[1]

"""

clusters = [0 1 1 ..., 1 0 0]
label[i] is the code or index of the centroid the i’th observation is closest to.

"""

cluster_with_protein_name = retrieve_protein(table, clusters_with_nb)
print cluster_with_protein_name