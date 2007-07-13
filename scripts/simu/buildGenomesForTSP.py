#! /users/ldog/muffato/python -OO

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myPhylTree

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("outFiles",str,""), ("geneLossRate",int,10), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
for anc in phylTree.listAncestr:
	genome = utils.myGenomes.EnsemblGenome(options["genesFile"] % phylTree.fileName[anc])
	ancGenes = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[anc])
	f = utils.myTools.myOpenFile(options["outFiles"] % phylTree.fileName[anc], "w")
	for c in genome.lstChr:
		random.shuffle(genome.lstGenes[c])
		genome.lstGenes[c] = genome.lstGenes[c][:(1.-options["geneLossRate"]/100)*len(genome.lstGenes[c])]
	for g in genome:
		(c,i) = ancGenes.dicGenes[g.names[0]]
		print >> f, g.chromosome, " ".join(ancGenes.lstGenes[c][i].names)
	f.close()

