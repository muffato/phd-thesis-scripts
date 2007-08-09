#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes:
		- Dessine la matrice des genes orthologues entre deux genomes.
		- Dessine le karyotype d'un genome face a l'autre
		- Renvoie les couples de chromosomes orthologues avec des stats sur l'evolution des genomes (nb genes, rearrangements)
		- Renvoie la liste des couples de genes orthologues avec les details sur leurs positions
		- Reordonne le genome 1 pour qu'il soit plus ressemblant au genome 2
		- Renvoie les diagonales entre les deux genomes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags
import utils.myPsOutput


#############
# FONCTIONS #
#############

# Fabrique la liste des orthologues entre les deux genomes
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
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
			res[c1][genome1.dicGenes[s][1]] = [(c,i) for (c,i) in tmp if c in chr2]
			
	nb12 = 0.
	nbTot = 0.
	synt = 0.
	for c1 in chr1:
		score = dict.fromkeys(chr2, 0)
		for i1 in res[c1]:
			for (c2,i2) in res[c1][i1]:
				score[c2] += 1
		for c2 in chr2:
			nb12 += (score[c2]*(score[c2]-1))/2.
		nbTot += (len(res[c1])*(len(res[c1])-1))/2.
		if len(res[c1]) != 0:
			synt += max(score.values())/float(len(res[c1]))

	print "%.4f\t%.4f" % (100.*nb12/nbTot, 100.*synt/len(chr1)),



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["referenceGenome"], [], __doc__)


genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])

for l in sys.stdin:
	# Chargement des fichiers
	genome1 = utils.myGenomes.loadGenome(l[:-1])
	print l[:-1],
	buildOrthosTable(genome1, genome1.lstChr, genome2, genome2.lstChr)
	buildOrthosTable(genome2, genome2.lstChr, genome1, genome1.lstChr)
	print

