
# -*- encoding: utf-8 -*-
# ** PROJET DATA MINING
# ** L.Leaute A.Souvane K.Kastano J.Estebeteguy


######################################################
# PACKAGES                                           #
######################################################

import sys,numpy as np,pandas as pd
from operator import itemgetter
from scipy.cluster.vq import vq, kmeans2, whiten
from scipy.cluster.hierarchy import linkage , fcluster, dendrogram
import matplotlib.pyplot as plt
import json

######################################################
# VARIABLES GLOBALES                                 #
######################################################

if len(sys.argv) >= 3:
	name = sys.argv[1]
	k = int(sys.argv[2]) # nombre de clusters kmeans
else:
	print "Usage: python script.py [fichier.tab] [nb] [1]"
	sys.exit()

######################################################
# PARAMETRES LES PLUS REPRESENTES                    #
###################################################### 

def representativity():
	''' 
	Calcule le pourcentage de presence d'un parametre pour l'ensemble des proteines 
	Renvoie un dictionnaire avec les parametres comme mot cles et leur pourcentages comme valeurs
	ex. {'Entry': 100.0, 'Length': 100.0, 'Mass': 100.0, 'GoBiologicalProcess': ['...'], 'GoCellularComponent': ['...'], 'GOMolecularFunction': ['...']}
	'''
		
	fic = open(name,"r")
	headers = fic.readline().rstrip().split('\t')
	dico = {}
	cptline = 0.0
	   
	for entry in headers:
		dico[entry]=0
	for line in fic.readlines():
		if not line :
		    break
		cptline +=1.0
		prot = line.rstrip().split('\t')
		for nb_col in range(len(prot)):
		    if prot[nb_col] != '':
		        dico[ headers[nb_col] ] += 1.0
		        
	for i in dico.keys():
		dico[i]= (dico[i]/cptline) * 100
		
	fic.close()
	
	return dico,headers,cptline



######################################################
# MISE EN FORME DU JEU DE DONNEE                     #
###################################################### 

def get_data():

	'''On implemente un dictionnaire qui contient une entree pour chaque protéine.
	Chaque protéine est un dictionnaire contenant les valeurs de chaque parametre.
	ex.
	{ 
  		PROT1 : { 'Length': '397', 'Mass': '45,318', 'GoBiologicalProcess':
            ['fructose 2,6-bisphosphate metabolic process
            [GO:0006003]'], 'GoCellularComponent': ['cytoplasm
            [GO:0005737]'], 'GOMolecularFunction':
            ['6-phosphofructo-2-kinase activity [GO:0003873]', 'ATP
            binding [GO:0005524]'] ,
  		PROT2 : { ... }
	}
	'''
	
	fic = open(name,"r")
	headers = fic.readline().rstrip().split('\t')
	
	data = {}
	go_cell = [] 
	go_mol = []
	go_bio = []

	for line in fic.readlines():
		if not line:
		    break
		prot = line.rstrip().split('\t')
		data[prot[0]] = {}

		for nb_col in range(1,len(prot)):

			# Mass
			if headers[nb_col] == "Mass":
				mass = prot[nb_col]
				mass = mass.replace(",","")
				data[prot[0]]["Mass"] = float(mass)
                        # Length
			elif headers[nb_col] == "Length":
				length = prot[nb_col]
				length = length.replace(",","")
				data[prot[0]]["Length"] = float(length)
			
			# Gene ontology (cellular component)
		        elif headers[nb_col] == "Gene ontology (cellular component)":
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


	# suppression des proteines sans go_terms

	for key in data.keys():
		if len(data[key]) != len(headers)-1:
			del data[key]
	
	return data,headers,go_cell,go_mol,go_bio

######################################################
# KMEANS SUR LA LONGUEUR ET LA MASSE                               
###################################################### 

def create_mass_length_array():
	'''
	Création d'un array qui contient les coordonnées de chaque protéine	
	'''

	list_tmp = []

	for prot in data:	    
		x = data[prot]['Mass']
		y = data[prot]['Length']
		list_tmp.append( [ x ,  y ] )

	return np.array( list_tmp )
	

def launch_kmeans(array_mass_length,k):

	whitened = whiten(array_mass_length) # normalisation des donnees division par l'ecart type en colonne
	centroids,cluster_index = kmeans2(whitened, k, iter=100,minit='points')

	return centroids,cluster_index

def create_2D_table(data):
	'''
	Ce tableau permettra de récuperer la protéine associée à chaque point de la méthode kmeans
	'''

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


def retrieve_protein(table, clusters_with_nb,k):
	clusters_with_name = [] 

	for i in range(k):
		clusters_with_name.append([])

	for i in range(len(clusters_with_nb)):
		protein_name = table[i][0]
		clusters_with_name[ clusters_with_nb[i] ].append(protein_name)

	return clusters_with_name



def retrieve_properties(un_cluster, go_category):
	"""
	Récupere le terme GO de la catégorie donnée qui est en commun pour toutes les protéines
	presentes dans un_cluster
	"""
	if un_cluster == []:
		return []

	goProt = []
	for prot in un_cluster:
		goProt.append(data[prot][go_category])

	result = set(goProt[0]).intersection(*goProt[:1])

	return list(result)



def statistics(cluster):
	nb_prot = len(cluster.proteins)
	mass_sum = 0
	length_sum = 0
	n = 0

	for prot in cluster.proteins:

		mass_sum = mass_sum + data[prot]['Mass']
		length_sum = length_sum + data[prot]['Length']
		n = n + 1

	mass_mean = mass_sum / n
	length_mean = length_sum / n

	return nb_prot, mass_mean, length_mean


######################################################
# MATRICE DES go_terms                                #
###################################################### 

def matrice_goterm(dico,go_type,go_list):
	'''
	Matrice binaire contenant pour chaque proteine (ligne) la presence ou l'absence de chaque goterm (colonne)
	'''
	matrice = []
	for prot in dico:
		row = []
		for goterm in go_list:
		    if goterm in dico[prot][go_type]:
		        row.append(1)
		    else:
		        row.append(0)        
		matrice.append(row)

	return matrice


######################################################
# REPRESENTATIVITE  DES go_terms                      #
###################################################### 

def go_representativite(dico,go_type,go_list):
	'''
	Calcul du pourcentage de representativite de chaque GoTerm pour voir si on peut supprimer certaines colonnes peu informatives.
	Le calcul de la matrice de dissimilarite sera moins long.
	'''
	gTerms = {}

	for go in go_list:
		gTerms[go]=0.0
	for prot in dico:
		for go in dico[prot][go_type]:
		    if go != '':
		        gTerms[go] += 1
	for i in gTerms.keys():
		gTerms[i]= (gTerms[i]/len(dico)) * 100
	
	gt = gTerms.items()
	gt.sort(key=itemgetter(1),reverse=True)
	return gTerms,gt


        
######################################################
# Clustering hiérarchique avec une catégorie GO
######################################################
        
def cluster_on_go_term(clustered_data, go_category, go_terms, data):
    """
    data: la structure de données initiale qui contient toutes les infos sur les protéines
    clustered_data: un cluster obtenu d'un clustering précedant, 
    c'est un tableau avec le noms de protéines qui font partie du cluster:
    ["prot1", "prot2", ..., "protn"]
    G0category: le nom de la catégorie GO utilisée pour le clustering
    go_terms: les termes récuperés du fichier de données initial
    
    Retourne: un tableau qui contient plusieurs tableaux qui correspondent à 
    tous les clusters des protéines produits:
    [["prot1,"prot2",...], ["prot3", "prot4", ...], ...]
    """

    if len(clustered_data) < 2:
        return []                
    
    scoresMat = []
    for protein in clustered_data:
        row =[]
        for goterm in go_terms:
            if goterm in data[protein][go_category]: 
                row.append(1.0)
            else:
                row.append(0.0)
        scoresMat.append(row)

    dist = linkage(scoresMat,'complete')
        
    # Tableau avec le numero du cluster auquel chaque protéine initiale appartient
    flatClusters = fcluster(dist,3,'distance')
    
    #cluster_output: un dataFrame avec une colonne protein qui contient toutes les protéines et 
    #une colonne cluster qui contient le numéro du cluster qui correspond à chacune de protéines
    cluster_output = pd.DataFrame({'prot':clustered_data, 'cluster':flatClusters})

    # Récuper tous les clusters de protéines de cluster_output
    clusterNames = []
    for cl in flatClusters:
        if cl not in clusterNames:
            clusterNames.append(cl)

    finalClusters = []
    for cl in clusterNames:
        clProt = cluster_output[cluster_output.cluster == cl]["prot"] .values
        finalClusters.append(clProt.tolist())

    return finalClusters

##########################################################################
# Structure d'abre pour stocker les résultats du clustering hiérarchique
##########################################################################

class tree_node:
	def __init__(self, _name, _proteins):
		self.name = _name
		self.proteins = _proteins
		self.children = []
		self.go_bio = retrieve_properties(self.proteins, "GoBiologicalProcess")
		self.go_cell = retrieve_properties(self.proteins, "GoCellularComponent")
		self.go_mol = retrieve_properties(self.proteins, "GoMolecularFunction")

	def print_cluster(self):
		print "--------------------------"
		print "Cluster", self.name, "\n"
		print "Proteines :", self.proteins,"\n"
		print "Termes GO Biological Process :", self.go_bio,"\n"
		print "Termes GO Cellular Component :", self.go_cell,"\n"
		print "Termes GO Molecular Function :", self.go_mol
		print "--------------------------"	


##############################################################
# Création de fichiers pour la visualisation sur Cytoscape   #
##############################################################

def stringify(table):
	res = ""
	for element in table:
		res += element+ " "

	return res

def write_cluster_info_in_file(fic, fic_info, cluster, cluster_parent):
	"""
	Écrit les informations sur le cluster donné dans les deux fichiers
	qui seront à utiliser sur Cytoscape
	"""

	fic.write(cluster_parent.name+"\tchild\t"+cluster.name+"\n")

	prot_string = stringify(cluster.proteins)
	fic_info.write(cluster.name+"\t"+prot_string+"\t")
		
	go_cell_string = stringify(cluster.go_cell)
	fic_info.write(go_cell_string+"\t")

	go_mol_string = stringify(cluster.go_mol)
	fic_info.write(go_mol_string+"\t")

	go_bio_string = stringify(cluster.go_bio)
	fic_info.write(go_bio_string+"\t")

	nb_prot, mass_mean, length_mean = statistics(cluster)
	fic_info.write(str(nb_prot)+"\t"+str(mass_mean)+"\t"+str(length_mean))

	fic_info.write("\n")

def create_cytoscape_files(root):

	fic = open("results.tab","w")
	fic_info = open("extra_info.tab", "w")
	fic.write("cluster1\tinteraction\tcluster2\n")
	fic_info.write("id\tproteins\tgo_cell\tgo_mol\tgo_bio\tnb_proteins\tmass_mean\tlength_mean\n")

	for i in root.children:

		fic.write("data\tchild\t"+i.name+"\n")
	
		for j in i.children:
			write_cluster_info_in_file(fic, fic_info, j, i)


			for k in j.children:
				write_cluster_info_in_file(fic, fic_info, k, j)
			

				for l in k.children:
					write_cluster_info_in_file(fic, fic_info, l, k)

	fic.close()
	fic_info.close()

############################################################
# Affichage alternatif
############################################################

def print_results(root):
	print "2"
	for i in root.children:
		i.print_cluster()
	
		for j in i.children:
			j.print_cluster()
		
			for k in j.children:
				k.print_cluster()

				for l in k.children:
					l.print_cluster()

######################################################
# MAIN                                               #
######################################################

# construction du jeu de donnees
data,headers,go_cell,go_mol,go_bio = get_data()


# Clustering Kmeans sur le poids et la taille 

array_mass_length = create_mass_length_array()
table = create_2D_table(data)
centroids,clusterIndex = launch_kmeans(array_mass_length,k)

kmeans_protein_clusters = retrieve_protein(table, clusterIndex,k)

# nombre de proteines par cluster kmeans
print "Clustering des kmeans sur le poids et la taille fini."

# contiendra les résultats finaux du clustering
clusters_tree = [] 

root = tree_node("data", data.keys())

print "Clustering sur les termes GO, patientez..."
for i in range(len(kmeans_protein_clusters)):


	lvl1 = cluster_on_go_term(kmeans_protein_clusters[i], "GoCellularComponent", go_cell, data)
	
	node = tree_node(str(i), kmeans_protein_clusters[i])
	root.children.append(node)

	clusters_tree.append([])
       
	for j in range(len(lvl1)):
	
		lvl2 = cluster_on_go_term(lvl1[j], "GoMolecularFunction" ,go_mol, data)
            
		index = str(i)+"."+str(j)
		node = tree_node(index, lvl1[j])
		root.children[i].children.append(node)	


		clusters_tree[i].append([])

               
		for k in range(len(lvl2)):
			lvl3 = cluster_on_go_term(lvl2[k], "GoBiologicalProcess",go_bio, data) 

			index = str(i)+"."+str(j)+"."+str(k)
			node = tree_node(index, lvl2[k])
			root.children[i].children[j].children.append(node)

			for l in range(len(lvl3)):
				index = str(i)+"."+str(j)+"."+str(k)+"."+str(l)
				node = tree_node(index, lvl3[l])
				root.children[i].children[j].children[k].children.append(node)
                        
			clusters_tree[i][j].append(lvl3)

print "Clustering fini."
create_cytoscape_files(root)
if (len(sys.argv)==4 and sys.argv[3] == "1"):
	print "1"
	print_results(root)

