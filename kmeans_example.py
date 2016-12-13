#coding:utf-8

from numpy import array
from scipy.cluster.vq import vq, kmeans2, whiten
import sys
import re

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
dico = {}
go_cell = [] 
go_mol = []
go_bio = []
k = 5 # number of cluster of kmeans method 

file_name = sys.argv[1]


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
				data[prot[0]]["GoMolecularFunction"] = gos
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
				terme = terme.replace("to ", "")
				terme = terme.replace("  ", " ")
				terme = re.sub(r"[*]", "", terme)
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

	termes_representatifs = [first, second, third]

	return termes_representatifs


def statistics(cluster):
	nb_prot = len(cluster)
	somme_mass = 0
	somme_length = 0
	n = 0

	for prot in cluster:
		somme_mass = somme_mass + float(data[prot]['Mass'].replace(',','')) 
		somme_length = somme_length + float(data[prot]['Length'].replace(',','')) 
		n = n + 1

	moyenne_mass = somme_mass / n
	moyenne_length = somme_length / n

	return nb_prot, moyenne_mass, moyenne_length


def creer_arbre_vide(racine, cluster, ett = None):
	sommets = [racine]
	fils = {racine: []}
	pere = {racine: None}
	etiquettes = { racine: ett }
	liste_proteine = {racine: cluster}
	return [ sommets, fils, pere, racine, etiquettes, liste_proteine ]

def ajouter_sommet(A,s,pere_s,proteines,etiquette=None):
	sommets,fils,pere,racine, ett, liste_proteine = A
	if s in sommets:
		raise ValueError ("Le sommet existe deja")
	if not pere_s in sommets:
		raise ValueError("Le pere n'existe pas dans l'arbre")

	sommets.append(s)
	fils[s] = []
	fils[pere_s].append(s)
	pere[s] = pere_s
	ett[s] = etiquette
	liste_proteine[s] = proteines

def etiquette(A,f):
	return A[4][f]

def racine(A):
	return A[3]

def sommets(A):
	return A[0]

def fils(A,s):
	return A[1][s]

def pere(A,s):
	if racine(A) == s:
		return None

	return A[2][s]

def proteines(A,s):
	return A[4][s]


def afficher_arbre(A):

	# Ouverture fichier
	fichier = open("visualisation_arbre.dot", "w") 
	
	# Ecriture dans fichier
	fichier.write("digraph G { graph [ ordering = ""out"" ];")

	# Ecriture Racine
	sommet = racine(A)
	ett = etiquette(A,racine(A))
	input0 = '{} [ label= "{}" ];'.format(sommet,ett)
	fichier.write(input0)

	# Ecriture tous les sommets
	afficher_arbre_recc(A,racine(A),fichier)
	fichier.write("}")

	# Fermeture fichier
	fichier.close()


def afficher_arbre_recc(A,s,fichier):

	for f in fils(A,s):
		sommet = f
		ett = etiquette(A,f)
		input1 = '{} [ label= "{}" ];'.format(sommet,ett)
		fichier.write(input1)

		input2 = '{} -> {};'.format(pere(A,sommet),sommet)
		fichier.write(input2)
		afficher_arbre_recc(A,f,fichier)


def contient_string(cluster):
    try:
    	for c in cluster:
    		re.search(r"P*", c)
        return True
    except TypeError:
        return False


def construire_arbre(cluster):
	A = creer_arbre_vide(0,cluster,"Cluster 0")
	construire_arbre_recc(A, cluster, 0, 0)
	return A

def construire_arbre_recc(A, cluster, sommet, pere):
	if contient_string(cluster) == True:
		return None

	for c in cluster:
		s = 0
		while s in sommets(A):
			s = s + 1 

		
		label = "Cluster " + str(s)
		ajouter_sommet(A, s, pere, c, label)
		sommet = s
		construire_arbre_recc(A, c, s, sommet)

	return A

def cluster_to_dico(cluster, cluster_nb):
	
	nb_prot, moyenne_mass, moyenne_length = statistics(cluster)
	process_representatifs = retrieve_properties(cluster, "GoBiologicalProcess")
	component_representatifs = retrieve_properties(cluster, "GoCellularComponent")
	function_representatifs = retrieve_properties(cluster, "GoMolecularFunction")

	dico[cluster_nb] = [ cluster, nb_prot, moyenne_mass, moyenne_length, process_representatifs, component_representatifs, function_representatifs ]


clusters = [	
				[
					['P51561', 'p474418'], ['P84845', 'P84845']
				], 
				
				[
					['P84845'], 
					['P84845', 'P84845', 'P84845', 'P84845']
				]
			]



create_data()

"""
data = 

{ 
  
  PROT1 : { 'Length': 397,
            'Mass': 45.318,
            'GoBiologicalProcess': ['fructose 2,6-bisphosphate metabolic process [GO:0006003]', 'fructose metabolic process [GO:0006000]', 'regulation of glycolytic process [GO:0006110]'], 
            'GoCellularComponent': ['cytoplasm [GO:0005737]'], 
            'GoMolecularFunction': ['6-phosphofructo-2-kinase activity [GO:0003873]', 'ATP binding [GO:0005524]']   ,

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

"""
for i in range(1,10):
	k = i
	kmean = launch_kmeans(array_mass_length)

	centroids = kmean[0]
	for j in range(len(centroids)):
		if centroids[j][0]<0 or centroids[j][1]<0:
			print "ok"
"""


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


clusters_with_protein_name = retrieve_protein(table, clusters_with_nb)


for i in range(len(clusters_with_protein_name)):
	print "CLUSTER ", i
	nb_prot, moyenne_mass, moyenne_length = statistics(clusters_with_protein_name[i])
	print "Nombre de protéines dans le cluster: ", nb_prot
	print "Moyenne des masses: ", moyenne_mass
	print "Moyenne des longueurs: ", moyenne_length
	print "Bio process: ", str(retrieve_properties(clusters_with_protein_name[i], "GoBiologicalProcess")).strip('[]')
	print "Cell comp: ", str(retrieve_properties(clusters_with_protein_name[i], "GoCellularComponent")).strip('[]')
	print "Mol function: ", str(retrieve_properties(clusters_with_protein_name[i], "GoMolecularFunction")).strip('[]')
	print "\n\n"


A = construire_arbre(clusters_with_protein_name)
afficher_arbre(A)

