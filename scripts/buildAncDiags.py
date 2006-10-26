#! /usr/bin/python2.4

__doc__ = """
Etant donne un ancetre dans l'arbre phylogenetique:
  Extrait toutes les diagonales de genes entre deux especes filles de branches differentes
  et d'une espece fille et d'une espece outgroup
  On reconstruit donc les ensembles de genes presents cote a cote chez l'ancetre.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("ancetre",str,""), ("output",str,""), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)
genesAncRoot = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.root, False)
genesAncNode = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % options["ancetre"], False)

# Les listes des especes entre lesquelles on cherche des diagonales
groupes = []
fils = set([])
for (e,_) in phylTree.items[options["ancetre"]]:
	tmp = phylTree.getSpecies(e)
	groupes.append(tmp)
	fils.update(tmp)
outgroup = set(listEspeces).difference(fils)

# Traduction des genomes en liste des genes ancestraux
genomesRoot = {}
for e in geneBank.dicEspeces:
	print >> sys.stderr, "Traduction de %s avec les genes %s ..." % (e, phylTree.root),
	genomesRoot[e] = utils.myDiags.translateGenome(geneBank.dicEspeces[e], genesAncRoot)
	print >> sys.stderr, "OK"
genomesNode= {}
for e in fils:
	print >> sys.stderr, "Traduction de %s avec les genes %s ..." % (e, options["ancetre"]),
	genomesNode[e] = utils.myDiags.translateGenome(geneBank.dicEspeces[e], genesAncNode)
	print >> sys.stderr, "OK"

# Liste des positions des genes ancestraux dans les genomes modernes
print >> sys.stderr, "Extraction des positions des genes de %s ..." % phylTree.root,
locationsRoot = utils.myDiags.buildAncGenesLocations(geneBank, genesAncRoot)
print >> sys.stderr, "de %s ..." % options["ancetre"],
locationsNode = utils.myDiags.buildAncGenesLocations(geneBank, genesAncNode)
print >> sys.stderr, "OK"


# La fonction qui permet de traiter les diagonales
def combinDiag(c1, c2, d1, d2):
	global combin, options
	global e1, e2, fils

#	if len(d1) != len(d2):
#		print >> sys.stderr, "PROBLEME.L"
#	if e1 in fils:
#		dodo1 = [genomesNode[e1][c1][i] for i in d1]
#	dada1 = [genomesRoot[e1][c1][i] for i in d1]
#	didi1 = [geneBank.dicEspeces[e1].lstGenes[c1][i].names[0] for i in d1]
#	if e2 in fils:
#		dodo2 = [genomesNode[e2][c2][i] for i in d2]
#	dada2 = [genomesRoot[e2][c2][i] for i in d2]
#	didi2 = [geneBank.dicEspeces[e2].lstGenes[c2][i].names[0] for i in d2]
#	if e1 in fils and e2 in fils:
#		if dodo1 != dodo2:
#			print >> sys.stderr, "PROBLEME.O"
#			sys.exit(0)
#	if dada1 != dada2:
#		print >> sys.stderr, "PROBLEME.A"
#		sys.exit(0)
	
	

	
	# Si on a demande les diagonales projetees sur un genome particulier
	if options["output"] == "":
		tmp = set([])
		if e1 in fils:
			tmp.update([genomesNode[e1][c1][i] for i in d1])
		if e2 in fils:
			tmp.update([genomesNode[e2][c2][i] for i in d2])
		d = [i for i in tmp if i != -1]
	else:
		if e1 == options["output"]:
			d = [geneBank.dicEspeces[e1].lstGenes[c1][i].names[0] for i in d1]
		elif e2 == options["output"]:
			d = [geneBank.dicEspeces[e2].lstGenes[c2][i].names[0] for i in d2]
		else:
			d = []

	if len(d) >= options["minimalLength"]:
		combin.addLink(d)

# On genere toutes les diagonales entre des paires d'especes qui passent par le noeud
combin = utils.myTools.myCombinator([])

for (i,j) in utils.myTools.myMatrixIterator(len(groupes), len(groupes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	for e1 in groupes[i]:
		for e2 in groupes[j]:
			print >> sys.stderr, "Extraction des diagonales (sous-especes) entre %s et %s ..." % (e1,e2),
			utils.myDiags.iterateDiags(genomesNode[e1], locationsNode[e2], options["fusionThreshold"], combinDiag)
			print >> sys.stderr, "OK"

for g in groupes:
	for e1 in g:
		for e2 in outgroup:
			print >> sys.stderr, "Extraction des diagonales (outgroup) entre %s et %s ..." % (e1,e2),
			utils.myDiags.iterateDiags(genomesRoot[e1], locationsRoot[e2], options["fusionThreshold"], combinDiag)
			print >> sys.stderr, "OK"

for g in combin.getGrp():
	print " ".join([str(x) for x in g])

