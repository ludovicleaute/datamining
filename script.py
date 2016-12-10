# ** PROJET DATA MINING
# ** L.Leaute A.Souvane K.Kastano J.Estebeteguy

######################################################
# PACKAGES                                           #
######################################################

import sys,numpy as np,pandas as pd
from operator import itemgetter
from scipy.cluster.vq import vq, kmeans2, whiten
from scipy.cluster.hierarchy import linkage , fcluster, dendrogram, cut_tree
import matplotlib.pyplot as plt

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
''' 
On implemente un dictionnaire qui contient en entrees toutes les proteines.
chaque prot est un dictionnaire contenant les valeurs de chaque parametre

{ 
  
  PROT1 : { 'Length': '397',
            'Mass': '45,318',
            'GoBiologicalProcess': ['fructose 2,6-bisphosphate metabolic process [GO:0006003]'],
            'GoCellularComponent': ['cytoplasm [GO:0005737]'], 
            'GOMolecularFunction': ['6-phosphofructo-2-kinase activity [GO:0003873]', 'ATP binding [GO:0005524]']   ,

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
  

# clustering hierarchique en fonction des compartiments


scoresMat = []
for protein in cluster_with_protein_name[0]:
	row =[]
   	for goterm in go_cell:
   		if goterm in data[protein]["GoCellularComponent"]: 
			row.append(1.0)
		else:
			row.append(0.0)
	scoresMat.append(row)

print "\nmatrice created successfully\n"

dist = linkage(scoresMat,'complete')

assignement = fcluster(dist,3,'distance')



cluster_output = pd.DataFrame({'prot':cluster_with_protein_name[0] , 'cluster':assignement})

print cluster_output[:20]

res = {}
for n in assignement:
	res[n]=0
for n in assignement:
	res[n]+=1

print [(k,v) for k,v in res.items() if v == max([v for v in res.values()])] 

print max([v for v in res.values()])
