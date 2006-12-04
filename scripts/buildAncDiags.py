#! /users/ldog/muffato/python

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
import utils.myMaths
import utils.myDiags


#############
# FONCTIONS #
#############

def calcDiags(e1, e2):

	# La fonction qui permet de traiter les diagonales
	def combinDiag(c1, c2, d1, d2):
		global diagEntry, geneBank

		if len(d1) < options["minimalLength"]:
			return
		
		dn1 = [geneBank.dicEspeces[e1].lstGenes[c1][i].names[0] for i in d1]
		dn2 = [geneBank.dicEspeces[e2].lstGenes[c2][i].names[0] for i in d2]

		for tmp in toStudy:
			diagEntry[tmp].append( ((e1,c1,dn1), (e2,c2,dn2)) )
			
	
	toStudy = [phylTree.getFirstParent(e1,e2)]
	for tmp in phylTree.items:
		s = phylTree.getSpecies(tmp)
		if (e1 in s) and (e2 not in s):
			toStudy.append( tmp )
		elif (e2 in s) and (e1 not in s):
			toStudy.append( tmp )

	# Chargement des orthologues
	genesAnc = utils.myGenomes.GenomeFromOrthosList(options["orthosFile"] % (e1,e2))
	nbGenesAnc = len(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr])
	del genesAnc.lstGenes

	newLoc = [[] for x in xrange(nbGenesAnc)]
	for e in [e1,e2]:
		genome = geneBank.dicEspeces[e]
		newGenome = {}
		for c in genome.lstChr:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
			if not options["keepOrthosLess"]:
				newGenome[c] = [x for x in newGenome[c] if x[0] != -1]
			if e == e2:
				for i in xrange(len(newGenome[c])):
					(ianc,s) = newGenome[c][i]
					if ianc != -1:
						newLoc[ianc].append( (c,i,s) )
		if e == e1:
			newGen = newGenome
	
	utils.myDiags.iterateDiags(newGen, newLoc, options["fusionThreshold"], options["sameStrand"], combinDiag)


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOrthosLess",bool,True), \
	("orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2")], \
	__doc__ \
)


# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)

# Pour sauver de la memoire
for esp in geneBank.dicEspeces:
	del geneBank.dicEspeces[esp].dicGenes

# La structure qui accueillera les diagonales et le calcul
diagEntry = dict( [(anc, []) for anc in phylTree.items] )
for (i,j) in utils.myTools.myMatrixIterator(len(listEspeces), len(listEspeces), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	print >> sys.stderr, '+',
	calcDiags(listEspeces[i], listEspeces[j])

# Plus besoin de ca ...
del geneBank

for anc in diagEntry:
	
	print >> sys.stderr, "Traitement de %s ..." % anc,
	lst = diagEntry[anc]
	print >> sys.stderr, len(lst),

	s = 0
	for ((e1,c1,d1),(e2,c2,d2)) in lst:
		s += len(d1)
		print "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (anc, len(d1), e1,c1," ".join(d1), e2,c2," ".join(d2))
	
	print >> sys.stderr, s, "OK"

