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

# Fabrique la liste des orthologues utilises pour le dessin
def buildOrthosTable(genome1, chr1, genome2, chr2):

	# Tous les orthologues entre pour les chromosomes OK du genome 1
	res = {}
	for c1 in chr1:
		res[c1] = {}
		for g1 in genome1.lstGenes[c1]:
			tmp = set()
			for s in g1.names:
				if s in genome2.dicGenes:
					tmp.add(genome2.dicGenes[s])
				if options["orthologuesList"] == "":
					continue
				if s in genesAnc.dicGenes:
					(c,i) = genesAnc.dicGenes[s]
					for ss in genesAnc.lstGenes[c][i].names:
						if ss in genome2.dicGenes:
							tmp.add(genome2.dicGenes[ss])
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
			# +/- includeGaps
			if (not options["includeGaps"]) and (len(tmp) == 0):
				continue
			res[c1][genome1.dicGenes[s][1]] = [(c,i) for (c,i) in tmp if c in chr2]
			

	return res


########
# MAIN #
########

# Arguments
modes = ["Matrix", "Karyotype", "OrthosChr"]
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps", bool, False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("output",str,modes), ("scaleY",bool,False), \
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("minHomology",int,90)], \
	__doc__
)


# Chargement des fichiers
#genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome1 = utils.myGenomes.AncestralGenome(noms_fichiers["studiedGenome"])
#genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
genome2 = utils.myGenomes.AncestralGenome(noms_fichiers["referenceGenome"])

nouveaux = 0
for g in genome1:


	for s in g.names:
		if s in genome2.dicGenes:
			(_,i) = genome2.dicGenes[s]
			print " ".join(genome2.lstGenes[utils.myGenomes.Genome.defaultChr][i].names)
			break
	else:
		print
		nouveaux += 1

sys.exit(0)
dupliques = 0
pertes = 0
for g in genome2:

	ens = set()
	for s in g.names:
		if s in genome1.dicGenes:
			ens.add(genome1.dicGenes[s])
	if len(ens) == 0:
		pertes += 1
	elif len(ens) > 1:
		print g.names
		dupliques += 1


print "Evolution de 1 vers 2"
print "nouveaux", nouveaux
print "dupliques", dupliques
print "pertes", pertes




