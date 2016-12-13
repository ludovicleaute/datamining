# -*- encoding: utf-8 -*-
# ** PROJET DATA MINING
# ** L.Leaute A.Souvane K.Kastano J.Estebeteguy


######################################################
# PACKAGES                                           #
######################################################

import sys,numpy as np,pandas as pd
from operator import itemgetter
from scipy.cluster.vq import vq, kmeans2, whiten
from scipy.cluster.hierarchy import linkage , fcluster, dendrogram#, cut_tree
import matplotlib.pyplot as plt
import json
import re


######################################################
# VARIABLES GLOBALES                                 #
######################################################

name = sys.argv[1]
k = int(sys.argv[2]) # nombre de clusters kmeans

######################################################
# PARAMETRES LES PLUS REPRESENTES                    #
###################################################### 
''' creation d'un dictionnaire qui contient les headers comme mot cles.
calcule du pourcentage d'information apporte par un parametre sur l'ensemble des proteines '''

def representativity():
		
	fic = open(name,"r") # ouverture du fichier pour recuperer les headers en global
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
	 
#>> {'Entry': 100.0, 'Length': 100.0, 'Mass': 100.0, 'GoBiologicalProcess': ['...'], 'GoCellularComponent': ['...'], 'GOMolecularFunction': ['...']}



######################################################
# MISE EN FORME DU JEU DE DONNEE                     #
###################################################### 
'''On implemente un dictionnaire qui contient en entrees toutes les proteines.
chaque prot est un dictionnaire contenant les valeurs de chaque parametre

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

def get_datas():
	
	fic = open(name,"r") # ouverture du fichier pour recuperer les headers en global
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


	# suppression des proteines sans goTerms

	for key in data.keys():
		if len(data[key]) != len(headers)-1:
			del data[key]
	
	return data,headers,go_cell,go_mol,go_bio

######################################################
# KMEANS LENGTH MASS                                 #
###################################################### 
'''
Mass_length_array is an array containing coordinates of each prot (as x = Mass and y = Length)
'''
def create_mass_length_array():

	list_tmp = []

	for prot in data:	    
		x = data[prot]['Mass']
		y = data[prot]['Length']
		list_tmp.append( [ x ,  y ] )

	return np.array( list_tmp )
	

def launch_kmeans(array_mass_length,k):

	whitened = whiten(array_mass_length) # normalisation des donnees division par l'ecart type en colonne
	centroids,clusterIndex = kmeans2(whitened, k, iter=100,minit='points')

	return centroids,clusterIndex

'''
This table will allow to retrieve the protein associated to each point of the kmeans method
'''
def create_2D_table(data):
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


def retrieve_properties(un_cluster, param):

    goProt = []
    for prot in un_cluster:
        goProt.append(data[prot][param])

    result = set(goProt[0]).intersection(*goProt[:1])

    return list(result)


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




######################################################
# MATRICE DES GOTERMS                                #
###################################################### 
'''
matrice binaire contenant pour chaque proteine (ligne) la presence ou l'absence de chaque goterm (colonne)
'''


def matrice_goterm(dico,go_type,go_list):
	matrice = []
	for prot in dico:
		row = []
		for goterm in go_list:
		    if goterm in dico[prot][go_type]:
		        row.append(1)
		    else:
		        row.append(0)        
		matrice.append(row)

	print "matrice created successfully\n"
	print "dimensions expected: ",len(dico),"/",len(go_list)
	print "dimensions obtained: ",len(matrice),"/",len(matrice[0]),"\n"

	return matrice


######################################################
# REPRESENTATIVITE  DES GOTERMS                      #
###################################################### 
'''
Calcul du pourcentage de representativite de chaque GoTerm pour voir si on peut supprimer certaines colonnes peu informatives
le calcul de la matrice de dissimilarite sera moins long
'''

def go_representativite(dico,go_type,go_list):
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
        
def clusterWithGOTerm(clustered_data, GOcategory, GOterms, data):
        """
        data: la structure de données initiale qui contient toutes les infos sur les protéines
        clustered_data: un cluster obtenu d'un clustering précedant, 
        c'est un tableau avec le noms de protéines qui font partie du cluster:
        ["prot1", "prot2", ..., "protn"]
        G0category: le nom de la catégorie GO utilisée
        GOterms: les termes récuperés du fichier de données initial
        
        Retourne: un tableau qui contient plusieurs tableaux qui correspondent à 
        tous les clusters des protéines produits:
        [["prot1,"prot2",...], ["prot3", "prot4", ...], ...]
        """
        if len(clustered_data) < 2:
                return []                
        
        scoresMat = []
        for protein in clustered_data:
	        row =[]
   	        for goterm in GOterms:
   		        if goterm in data[protein][GOcategory]: 
			        row.append(1.0)
		        else:
			        row.append(0.0)
	        scoresMat.append(row)

        #print "\nmatrice created successfully\n"

        dist = linkage(scoresMat,'complete')

        # Tableau avec le numero du cluster auquel chaque protéine initiale appartient
        flatClusters = fcluster(dist,3,'distance')


        #cluster_output: un dataFrame avec une colonne protein qui contient toutes les protéines et une colonne cluster qui contient le numéro du cluster qui correspond à chacune de protéines
        cluster_output = pd.DataFrame({'prot':clustered_data, 'cluster':flatClusters})

        # Avoir tous les différents clusters
        clusterNames = []
        for cl in flatClusters:
                if cl not in clusterNames:
                        clusterNames.append(cl)

        # Selectionner les protéines pour chaque cluster
        finalClusters = []
        for cl in clusterNames:
                clProt = cluster_output[cluster_output.cluster == cl]["prot"].values
                finalClusters.append(clProt.tolist())

        #print finalClusters, "\n----------------------"
        return finalClusters




######################################################
# MAIN                                               #
######################################################

# construction du jeu de donnees
data,headers,go_cell,go_mol,go_bio = get_datas()

# clustering des kmeans sur le poids et la taille
array_mass_length = create_mass_length_array()
table = create_2D_table(data)
centroids,clusterIndex = launch_kmeans(array_mass_length,k)

cluster_with_protein_name = retrieve_protein(table, clusterIndex,k)

# nombre de proteines par cluster kmeans

for i in range(k):
	print "nb proteines cluster",i," : ",len(cluster_with_protein_name[i])
'''
# Un peu de visualisation 
dat = pd.DataFrame(whiten(array_mass_length))	
coord = dat.as_matrix(columns=[0,1])


plt.figure(figsize=(10, 10), dpi=100)
plt.scatter(coord[:,0], coord[:,1], c=clusterIndex, s=25, cmap="winter")
plt.scatter(centroids[:,0],centroids[:,1],c="r",s=50)
plt.show()

# representativite des goterms
gterms,gt = go_representativite(data,"GoCellularComponent",go_cell)

print "\nGo-Terms representativity :\n"

for k,v in gt[:5]:
    print k 
    print v,"%"

print "\n"
'''
  


"""
DIANA:
4 niveaux de clustering:
- Clustering avec l'algorithme kmeans (en haut) 
- Clustering avec le terme GO Cellular Component
- Clustering avec le terme GO Biological Process
- Clustering avec le terme GO Molecular Function
Résultat: des clusters dans des clusters dans des clusters...
"""

clustersTree = [] # contiendra les résultats finaux du clustering

for i in range(len(cluster_with_protein_name)):
        #on ajoute un tableau pour chaque cluster du premier niveau
        clustersTree.append([])
        
        lvl1 = clusterWithGOTerm(cluster_with_protein_name[i], "GoCellularComponent", go_cell, data)
        

        for j in range(len(lvl1)):
                #on ajoute un tableau dans le tableau du premier niveau pour chaque cluster du deuxième niveau
                clustersTree[i].append([])

                lvl2 = clusterWithGOTerm(lvl1[j], "GoMolecularFunction",go_bio, data)

                for k in range(len(lvl2)):
                        
                        lvl3 = clusterWithGOTerm(lvl2[k], "GoBiologicalProcess",go_mol, data)
        
                        clustersTree[i][j].append(lvl3)



# Visualisation d'une partie des résultats:
x = 0
dico_arbre = {}
for i in range(len(clustersTree)):
    #print "UN CLUSTER COMMENCE.."
    for j in range(len(clustersTree[i])):
        #print "Un cluster dans le cluster"
        for k in range(len(clustersTree[i][j])):
            for f in range(len(clustersTree[i][j][k])):
                #print "un cluster dans le cluster dans le cluster"
                #print clustersTree[i][j][k]
                nb = str(i) + " " + str(j) + " " + str(k) + " " + str(f)
                dico_arbre[nb] = clustersTree[i][j][k][f]
                x+= len(clustersTree[i][j][k][f])
                


""""
# Ne marche pas à cause de l'array issu de DataFrame avec .values :/
with open("clustersTree.json", "w") as outfile:
    json.dump(clustersTree, outfile)

print("file saved")

"""       

"""
for i in range(len(clustersTree)):
	print "CLUSTER ", i
	#nb_prot, moyenne_mass, moyenne_length = statistics(clustersTree[i])
	#print "Nombre de protéines dans le cluster: ", nb_prot
	#print "Moyenne des masses: ", moyenne_mass
	#print "Moyenne des longueurs: ", moyenne_length
	print "Bio process: ", str(retrieve_properties(clustersTree[i], "GoBiologicalProcess")).strip('[]')
	print "Cell comp: ", str(retrieve_properties(clustersTree[i], "GoCellularComponent")).strip('[]')
	print "Mol function: ", str(retrieve_properties(clustersTree[i], "GoMolecularFunction")).strip('[]')
	print "\n\n"
"""

for num_cluster, cluster in dico_arbre.items():
    print "CLUSTER ", num_cluster
    #nb_prot, moyenne_mass, moyenne_length = statistics(clustersTree[i])
    #print "Nombre de protéines dans le cluster: ", nb_prot
    #print "Moyenne des masses: ", moyenne_mass
    #print "Moyenne des longueurs: ", moyenne_length
    print "Bio process: ", str(retrieve_properties(cluster, "GoBiologicalProcess"))
    print "Cell comp: ", str(retrieve_properties(cluster, "GoCellularComponent"))
    print "Mol function: ", str(retrieve_properties(cluster, "GoMolecularFunction"))
    print "nombre de proteines: ",len(cluster)
    print "\n\n"
