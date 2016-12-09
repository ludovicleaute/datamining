#coding:utf-8

from numpy import array
from scipy.cluster.vq import vq, kmeans2, whiten

"""

- faire structure de données avec tout (done)
- faire tableau 2D avec nom prot et [taille+longueur] (done)
- faire une liste de cluster avec prot à l'intérieur (done)
- calculer moyenne et stat dessus, ça peut nous permettre de trouver des outliers.... z-score permet de dire si une donnée est utilisable ou non. A FAIRE ABSOLUMENT 



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


"""


data = {}
go_cell = [] 
go_mol = []
go_bio = []
k = 4 # number of cluster of kmeans method 

file_name = "data.csv"


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


			# Gene ontology (cellular component)
			if headers[nb_col] == "Gene ontology (cellular component)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GoCellularComponent"] = gos
				for go in gos:
					if go not in go_cell:
						go_cell.append(go)


			# Gene ontology (molecular function)
			elif headers[nb_col] == "Gene ontology (molecular function)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GoMolecularFunction"] = gos
				for go in gos:
					if go not in go_mol:
						go_mol.append(go)
 

			# Gene ontology (biological process)
			elif headers[nb_col] == "Gene ontology (biological process)":
				gos = prot[nb_col].split('; ')
				data[prot[0]]["GoBiologicalProcess"] = gos
				for go in gos:
					if go not in go_bio:
						go_bio.append(go)


			else:
				data[prot[0]][headers[nb_col]]=prot[nb_col]

	fic.close()



	"""
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
				gos = gos.replace(",","")
				data[prot[0]]["Mass"] = float(gos)


            # Length
			if headers[nb_col] == "Length":
				gos = prot[nb_col]
				gos = gos.replace(",","")
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
	"""

# Mass_length_array is an array containing coordinates of each prot (as x = Mass and y = Length)
def create_mass_length_array():

	list_tmp = []

	for prot in data:	    
		x = float(data[prot]['Length'].replace(",",""))
		y = float(data[prot]['Mass'].replace(",",""))
		list_tmp.append( [ x ,  y ] )

	return array( list_tmp )
	

def launch_kmeans(array_mass_length):

	# Whithened normalize data...je sais pas ce que ça veut dire mais il faut le faire
	whitened = whiten(array_mass_length)
	kmean = kmeans2(whitened,k, 100, 1e-05, 'points')

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



def retrieve_properties(un_cluster, param):
	first = ""
	second = ""
	third = ""

	dico = {}
	for prot in un_cluster:
		if param in data[prot].keys():
			terme_liste = data[prot][param]
			for terme in terme_liste:
				terme = terme.replace("of ", "")
				terme = terme.replace("from ", "")
				terme = terme.replace("  ", "")
				terme_tab = terme.split(" ")
				for mot in terme_tab:
					if mot not in dico.keys():
						dico[mot] = 1
					else: 
						dico[mot] = dico[mot] + 1
	max1 = 0
	max2 = 0
	max3 = 0

	for mot in dico.keys():
		if dico[mot] > max1:
			max1 = dico[mot]
			first = mot
		elif dico[mot] > max2:
			max2 = dico[mot]
			second = mot
		elif dico[mot] > max3:
			max3 = dico[mot]
			third = mot

	return first, second, third





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
for i in range(1,10):
	k = i
	kmean = launch_kmeans(array_mass_length)

	centroids = kmean[0]
	for j in range(len(centroids)):
		if centroids[j][0]<0 or centroids[j][1]<0:
			print "ok"

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
#print cluster_with_protein_name

print retrieve_properties(cluster_with_protein_name[0], "GoBiologicalProcess")
