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


# Renvoie la liste des chromosomes orthologues
def getOrthosChr(table, chr):
	
	res = {}
	for c1 in chr:

		# On ne garde que les chromosomes orthologues majoritaires
		lst = [x for (x,_) in utils.myMaths.flatten(table[c1].itervalues())]
		count = [(lst.count(x),x) for x in set(lst)]
		count.sort()
		nb = (len(lst)*options["minHomology"])/100
		tmp = []
		for (n,c2) in count.__reversed__():
			tmp.append( (c2,n) )
			nb -= n
			if nb <= 0:
				break
		res[c1] = tmp
	return res




########
# MAIN #
########

# Arguments
modeMatrix = "Matrix"
modeKaryo = "Karyotype"
modeGenomeEvol = "GenomeEvolution"
modeOrthos = "OrthosGenes"
modeReindexedChr = "ReindexedChr"
modeDiags = "Diags"
modes = [modeMatrix, modeKaryo, modeOrthos, modeGenomeEvol, modeReindexedChr, modeDiags]
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("output",str,modes), ("reverse",bool,False), ("scaleY",bool,False), ("minHomology",int,90), \
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.loadGenome(options["orthologuesList"])
else:
	genesAnc = genome2
try:
	colors = utils.myGenomes.loadGenome(options["colorFile"])
except Exception:
	colors = options["defaultColor"]

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)




nb12 = 0.
nbTot = 0.
synt = 0.
table12 = buildOrthosTable(genome1, chr1, genome2, chr2)
for c1 in chr1:
	score = dict.fromkeys(chr2, 0)
	for i1 in table12[c1]:
		for (c2,i2) in table12[c1][i1]:
			score[c2] += 1
	for c2 in chr2:
		nb12 += (score[c2]*(score[c2]-1))/2.
	nbTot += (len(table12[c1])*(len(table12[c1])-1))/2.
	synt += max(score.values())/float(len(table12[c1]))
print 100.*nb12/nbTot
print 100.*synt/len(chr1)
