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
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)


# Pour sauver de la memoire
for esp in geneBank.dicEspeces:
	del geneBank.dicEspeces[esp].dicGenes

# La structure qui accueillera les diagonales
diagEntry = dict( [(anc, utils.myTools.myCombinator([])) for anc in phylTree.items] )

# La fonction qui permet de traiter les diagonales
def combinDiag2(c1, c2, d1, d2):
	global diagEntry, options, toStudy
	global e1, e2, anc, genomes

	if len(d1) < options["minimalLength"]:
		return
	
	# On rajoute la diagonale a l'ancetre commun et a chaque ancetre intermediaire
	
	#dd1 = [genomes[anc][e1][c1][i] for i in d1]
	#dd2 = [genomes[anc][e2][c2][i] for i in d2]
	#if dd1 != dd2:
	#	print >> sys.stderr, "PB !!!", anc, e1, e2, d1, d2, dd1, dd2
	#diagEntry[anc].addLink(dd1)

	for (e, tmp) in toStudy:
		if e == e1:
			dd = [genomes[tmp][e][c1][i] for i in d1]
		else:
			dd = [genomes[tmp][e][c2][i] for i in d2]
		diagEntry[tmp].addLink(dd)
			

# 2. On prepare tous les genomes ancestraux, les genomes traduits ...
genomes = {}
locations = {}
for anc in phylTree.items:
	
	genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc, False, False)

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
	locations[anc] = utils.myDiags.buildAncGenesLocations(geneBank, genesAnc, fils)
	print >> sys.stderr, "OK"
	
	del genesAnc

del geneBank

for anc in phylTree.items:

	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	print >> sys.stderr, "Extraction des diagonales de %s " % anc,

	for (i,j) in utils.myTools.myMatrixIterator(len(groupes), len(groupes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		for e1 in groupes[i]:
			for e2 in groupes[j]:

				toStudy = [(e1,anc)]
				for tmp in phylTree.items:
					s = phylTree.getSpecies(tmp)
					if (e1 in s) and (e2 not in s):
						toStudy.append( (e1,tmp) )
					elif (e2 in s) and (e1 not in s):
						toStudy.append( (e2,tmp) )
			
				utils.myDiags.iterateDiags(genomes[anc][e1], locations[anc][e2], options["fusionThreshold"], combinDiag2)
				sys.stderr.write(".")
	del locations[anc]
	print >> sys.stderr, " OK"

del genomes

for anc in diagEntry:
	for g in diagEntry[anc]:
		print anc, " ".join([str(x) for x in g])

