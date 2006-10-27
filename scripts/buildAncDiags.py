#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales entre chaque paire d'especes.
Les diagonales apportent les genes qui etaient sur un meme chromosome
  depuis leur ancetre commun dans les deux lignees.
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
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)

# La structure qui accueillera les diagonales
diagEntry = dict( [(anc, utils.myTools.myCombinator([])) for anc in phylTree.items] )

# La fonction qui permet de traiter les diagonales
def combinDiag2(c1, c2, d1, d2):
	global diagEntry, options
	global e1, e2, anc, genomes

	if len(d1) < options["minimalLength"]:
		return
	
	# On rajoute la diagonale a l'ancetre commun et a chaque ancetre intermediaire
	
	dd1 = [genomes[anc][e1][c1][i] for i in d1]
	dd2 = [genomes[anc][e2][c2][i] for i in d2]
	if dd1 != dd2:
		print >> sys.stderr, "PB !!!", anc, e1, e2, d1, d2, dd1, dd2
	diagEntry[anc].addLink(dd1)

	for (tmp,_) in phylTree.items:
		s = phylTree.getSpecies(tmp)
		if ((e1 in s) and (e2 not in s)) or ((e1 not in s) and (e2 in s)):
			dd1 = [genomes[tmp][e1][c1][i] for i in d1]
			dd2 = [genomes[tmp][e2][c2][i] for i in d2]
			if dd1 != dd2:
				print >> sys.stderr, "PB !!!", tmp, e1, e2, d1, d2, dd1, dd2
			diagEntry[tmp].addLink(dd1)
			

# 2. On prepare tous les genomes ancestraux, les genomes traduits ...
genomes = {}
locations = {}
for anc in phylTree.items:
	
	genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc, False)

	# Les listes des especes entre lesquelles on cherche des diagonales
	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	fils = utils.myMaths.flatten(groupes)
	
	# Traduction des genomes en liste des genes ancestraux
	print >> sys.stderr, "Traduction avec les genes de", anc, "",
	genomes[anc] = {}
	for e in fils:
		genomes[anc][e] = utils.myDiags.translateGenome(geneBank.dicEspeces[e], genesAnc)
		sys.stderr.write(".")
	print >> sys.stderr, " OK"
	
	# Liste des positions des genes ancestraux dans les genomes modernes
	print >> sys.stderr, "Extraction des positions des genes de %s ..." % anc,
	locations[anc] = utils.myDiags.buildAncGenesLocations(geneBank, genesAnc)
	print >> sys.stderr, "OK"
	
	del genesAnc



for anc in phylTree.items:

	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	print >> sys.stderr, "Extraction des diagonales de", anc,

	for (i,j) in utils.myTools.myMatrixIterator(len(groupes), len(groupes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		for e1 in groupes[i]:
			for e2 in groupes[j]:
				utils.myDiags.iterateDiags(genomes[anc][e1], locations[anc][e2], options["fusionThreshold"], combinDiag2)
				sys.stderr.write(".")
	print >> sys.stderr, " OK"



sys.exit(0)

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

