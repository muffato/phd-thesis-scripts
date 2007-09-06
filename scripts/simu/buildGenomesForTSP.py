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
(noms_fichiers, options) = utils.myTools.checkArgs( ["genesFile", "ancGenesFile"], [("geneLossRate",float,10)], "Melange un genome et enleve des genes aleatoirement" )

genome = utils.myGenomes.Genome(noms_fichiers["genesFile"])
for c in genome.lstChr:
	for i in xrange(10):
		random.shuffle(genome.lstGenes[c])
	genome.lstGenes[c] = genome.lstGenes[c][:int((1.-options["geneLossRate"]/100)*len(genome.lstGenes[c]))]

ancGenes = utils.myGenomes.Genome(noms_fichiers["ancGenesFile"])
for g in genome:
	(c,i) = ancGenes.dicGenes[g.names[0]]
	print g.chromosome, " ".join(ancGenes.lstGenes[c][i].names)

