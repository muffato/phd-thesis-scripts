#! /users/ldog/muffato/python -OO

__doc__ = """
Dessine la matrice des genes orthologues entre deux genomes.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPsOutput


#############
# FONCTIONS #
#############

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["referenceGenome", "studiedGenome"], [], __doc__ )


# Chargement des fichiers
genome2 = utils.myGenomes.AncestralGenome(noms_fichiers["studiedGenome"])
genome1 = utils.myGenomes.AncestralGenome(noms_fichiers["referenceGenome"])

nouveaux = 0
for g in genome2:

	for s in g.names:
		if s in genome1.dicGenes:
			(_,i) = genome1.dicGenes[s]
			#print " ".join(genome1.lstGenes[utils.myGenomes.Genome.defaultChr][i].names)
			break
	else:
		#print
		nouveaux += 1

dupliques = {}
pertes = 0
for g in genome1:

	ens = set()
	for s in g.names:
		if s in genome2.dicGenes:
			ens.add(genome2.dicGenes[s])
	if len(ens) == 0:
		pertes += 1
	elif len(ens) > 1:
		#print g.names
		dupliques[len(ens)] = dupliques.get(len(ens), 0) + 1


print "Evolution de A (%d genes) vers B (%d genes)" % \
	(len(genome1.lstGenes[utils.myGenomes.Genome.defaultChr]), len(genome2.lstGenes[utils.myGenomes.Genome.defaultChr]))
print "nouveaux", nouveaux
print "dupliques", sum(dupliques.values()), dupliques
print "pertes", pertes




