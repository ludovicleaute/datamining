# Projet de Datamining 

Ce projet est réalisé par Kristina Kastano, Julien Estebeteguy, Alexia Souvane et Ludovic Léauté</br>
L'objectif est de répondre à la question très générales :  quelles sont les protéines semblables ?</br>
Pour y parvenir nous avons travaillé sur des données issues de la bdd en ligne Swissprot.

Le projet est écrit en python 2.7.

# Instructions

<ul>
<li>Prenez bien soins d'installer les différents packages de python nécessaires et détaillés au début du fichier script.py</li>
<li>Pour réaliser le clustering executez dans un terminal la commande :</br>
> python script.py [fichier.tab] [nb] [1]</br>
fichier.tab est le fichier avec les protéines
Nb sera le nombre de clusters choisi pour le clustering de k-means.</br>
Le dernier chiffre (optionel) indique si vous voulez une visualization des résultats sur la console.
<li>Après l'éxecution vous pourrez utilisez les deux fichiers tab créés pour visualiser les résultats sur Cytoscape:
	<ul>
	<li> Importez le fichier results.tab (File > import > network > file). Choisez le cluster1 comme Source node et le cluster2 comme Target node</li>
	<li> Vous pouvez changer le layout en Organic (Layout > yFilesLayouts > Organic) pour une visualization plus claire</li>
	<li> Importez le fichier extra_info.tab (File > import > table) pour integrez les informations sur chaque noeud du réseau. </li>
	<li> La visualization est prête ! Vous pouvez cliquer sur un noeud pour voir les informations pertinantes affichées en bas.</li>
	</ul>
</ul>
