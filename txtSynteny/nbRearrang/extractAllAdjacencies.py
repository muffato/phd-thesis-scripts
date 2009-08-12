#! /users/ldog/muffato/python

__doc__ = """
	A partir d'un genome de reference et de 3 genomes voisins, extrait la liste des rearrangements specifiques de chaque branche
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("genomes",utils.myTools.FileList(2))], \
	[("minLength",int,1), ("eliminateDup",bool,False), ("removeNestedGenes",bool,False)], \
	__doc__ \
)

listEsp = arguments["genomes"]

# Genomes
genome = {}
for esp in listEsp:
	genome[esp] = utils.myGenomes.Genome(esp)

# Genes orthologues en commun sur des chromosomes
orthos = []
for gene in genome[listEsp[0]]:
	res = {}
	for esp in listEsp:
		pos = [(c,i) for (c,i) in genome[esp].getPosition(gene.names) if (c in genome[esp].chrSet[utils.myGenomes.ContigType.Chromosome]) or (c in genome[esp].chrSet[utils.myGenomes.ContigType.Scaffold])]
		if arguments["eliminateDup"] and (len(pos) != 1):
			break
		if len(pos) == 0:
			break
		res[esp] = pos
	else:
		orthos.append(res)

# Nombre de genes par chromosome
def doCount(esp):
	count = collections.defaultdict(int)
	for fam in orthos:
		for (c,_) in fam[esp]:
			count[c] += 1
	return count

# Nouvelle liste d'orthologues sur des scaffolds >= minLength
def rebuildOrthos():
	modif = False
	neworthos = []
	for fam in orthos:
		res = {}
		for esp in listEsp:
			pos = [(c,i) for (c,i) in fam[esp] if counts[esp][c] >= arguments["minLength"]]
			modif = modif or (len(pos) != len(fam[esp]))
			if len(pos) == 0:
				break
			res[esp] = pos
		else:
			neworthos.append(res)
	return (neworthos,modif)


# Boucle de filtrage des genes
modif = True
while modif:

	print >> sys.stderr, len(orthos), "genes"

	counts = {}
	for esp in listEsp:
		counts[esp] = doCount(esp)

	(orthos,modif) = rebuildOrthos()

# Nouveaux genomes restreints aux genes orthologues communs
def rewriteGenome(esp):
	lstGenes = collections.defaultdict(list)
	for (i,fam) in enumerate(orthos):
		for (c,x) in fam[esp]:
			s = genome[esp].lstGenes[c][x].strand
			lstGenes[c].append((x,i,s))
	adj = set()
	extr = set()
	for c in lstGenes:
		lstGenes[c] = [(i,s) for (_,i,s) in sorted(lstGenes[c])]
		assert len(lstGenes[c]) >= arguments["minLength"]
		for ((i1,s1),(i2,s2)) in utils.myTools.myIterator.slidingTuple(lstGenes[c]):
			if i1 < i2:
				adj.add( ((i1,s1),(i2,s2)) )
			elif i2 < i1:
				adj.add( ((i2,-s2),(i1,-s1)) )
		extr.add( lstGenes[c][-1] )
		(i1,s1) = lstGenes[c][0]
		extr.add( (i1,-s1) )
	return (lstGenes,adj,extr)

# Reduction des genomes
newGenome = {}
adjacencies = {}
extremities = {}
for esp in listEsp:
	(newGenome[esp],adjacencies[esp],extremities[esp]) = rewriteGenome(esp)
	print >> sys.stderr, "%s: %d chrom / %d adj / %d extr / %d genes" % (esp,len(newGenome[esp]),len(adjacencies[esp]),len(extremities[esp]),sum([len(x) for x in newGenome[esp].itervalues()]))

# Le jeu total d'adjacences
allAdj = set()
for esp in listEsp:
	allAdj.update(adjacencies[esp])

# Conserve ou casse
def status(esp, ((i1,s1), (i2,s2))):
	if ((i1,s1), (i2,s2)) in adjacencies[esp]:
		return "+"
	if ((i1,s1) in extremities[esp]) and ((i2,-s2) in extremities[esp]):
		return "?"
	return "-"

# Tableau final
for x in allAdj:
	def names(i, esp):
		return "/".join(genome[esp].lstGenes[c][x].names[0] for (c,x) in orthos[i][esp])
	r = [status(esp, x) for esp in listEsp] + ["%d/%d" % (x[0][1],x[1][1])] + [names(x[0][0], esp) for esp in listEsp]+ [names(x[1][0], esp) for esp in listEsp]
	print "\t".join(r)

