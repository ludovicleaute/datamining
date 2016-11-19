# ** PROJET DATA MINING
# ** L.Leaute A.Souvane K.Kastano J.Estebeteguy

######################################################
# PACKAGES                                           #
######################################################

import numpy as np
from operator import itemgetter
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

######################################################
# LECTURE DES DONNEES                                #
######################################################

name = "uniprot.tab" 
fic = open(name,"r")
headers = fic.readline().rstrip().split('\t') #noms des colonnes


######################################################
# PARAMETRES LES PLUS REPRESENTES                    #
###################################################### 
''' creation d'un dictionnaire qui contient les headers comme mot cles.
calcule du pourcentage d'information apporte par un parametre sur l'ensemble des proteines '''

'''
dico={}
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

print(dico)

#>> {'Entry': 100.0, 'Length': 100.0, 'Mass': 100.0, 'Gene ontology IDs': 88.46897783067996}

'''

######################################################
# MISE EN FORME DU JEU DE DONNEE                     #
###################################################### 
''' 
On implemente un dictionnaire qui contient en entrees toutes les proteines.
chaque prot est un dictionnaire contenant les valeurs de chaque parametre
{PROT1 : { 'Mass': '12', 'Length': '12', 'GoTerms' :[G0-0000001,GO-0000002] }, PROT2 :{...}; ...}
'''

data = {}
go_list = [] #liste de tous les goterms existants dans le jeu de donnees

for line in fic.readlines():
	if not line:
		break
	prot = line.rstrip().split('\t')
	data[prot[0]] = {}
	for nb_col in range(1,len(prot)):
		if headers[nb_col] == "Gene ontology IDs":
			gos = prot[nb_col].split('; ')
			data[prot[0]]["GoTerms"] = gos
			for go in gos:
				if go not in go_list:
					go_list.append(go)
		else:
			data[prot[0]][headers[nb_col]]=prot[nb_col]
        
fic.close()

# suppression des proteines sans goTerms

for key in data.keys():
	if len(data[key]) != len(headers)-1:
		del data[key]


######################################################
# MATRICE DES GOTERMS                                #
###################################################### 
'''
matrice binaire contenant pour chaque proteine (ligne) la presence ou l'absence de chaque goterm (colonne)
'''

dimCol = len(go_list)
dimRow = len(data)

'''
m = []

for prot in data:
    row = []
    for goterm in go_list:
        if goterm in data[prot]["GoTerms"]:
            row.append(1)
        else:
            row.append(0)        
    m.append(row)

print "matrice created successfully\n"

print "dimensions expected: ",dimRow,"/",dimCol
print "dimensions obtained: ",len(m),"/",len(m[0]),"\n"
'''

'''
Calcul du pourcentage de representativite de chaque GoTerm pour voir si on peut supprimer certaines colonnes peu informatives
le calcul de la matrice de dissimilarite sera moins long
'''

gTerms = {}

for go in go_list:
    gTerms[go]=0.0
for prot in data:
    for go in data[prot]["GoTerms"]:
        if go != '':
            gTerms[go] += 1

for i in gTerms.keys():
    gTerms[i]= (gTerms[i]/dimCol) * 100

gt = gTerms.items()
gt.sort(key=itemgetter(1),reverse=True)


print "Go-Terms representativity : \n"

for k,v in gt[:5]:
    print k 
    print v,"%"


'''
On peut choisir arbitrairement de ne garder que les 100 Goterms les plus representes cad les Goterms presents dans 60 proteines minimum (environs 1% de représentativité) le max etant 23%...
'''

bestGos = []

for k,v in gt[:100]:
    bestGos.append(k)
  
mbest = []    


for prot in data:
    row = []
    for goterm in bestGos:
        if goterm in data[prot]["GoTerms"]:
            row.append(1.0)
        else:
            row.append(0.0)        
    mbest.append(row)

print "\nmatrice created successfully\n"


print "dimensions expected: ",len(data),"/",len(bestGos)
print "dimensions obtained: ",len(mbest),"/",len(mbest[0]),"\n"

######################################################
# MATRICE DE DISTANCE ET CLUSTERING                  #
######################################################


clust = hierarchy.linkage(mbest,'average')

plt.figure()

dn = hierarchy.dendrogram(clust)

plt.show()

