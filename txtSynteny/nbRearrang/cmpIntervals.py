#!/usr/bin/env python2

__doc__ = """
	Parcourt le genome moderne et definit chacun des ses intervalles en terme de degre de conservation
"""

import sys
import itertools

import utils.myTools
import utils.myGenomes


# Argument:
arguments = utils.myTools.checkArgs( [("ancGenome",file), ("modernGenome",file)], [("minimalLength",int,0)], __doc__)

ancGenome = utils.myGenomes.Genome(arguments["ancGenome"])
genome = utils.myGenomes.Genome(arguments["modernGenome"])

# Reecriture des genomes
def rewriteGenome(genome):
	newGenome = {}
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		if len(genome.lstGenes[chrom]) >= abs(arguments["minimalLength"]):
			newGenome[chrom] = [(gene.names[0],gene.strand) for gene in genome.lstGenes[chrom]]
	return newGenome

newGenome = rewriteGenome(genome)
print >> sys.stderr, "modernGenome", sum([len(x) for x in newGenome.itervalues()])

newAncGenome = rewriteGenome(ancGenome)
print >> sys.stderr, "ancGenome", sum([len(x) for x in newAncGenome.itervalues()])

# Table de conversion entre les noms des deux genomes
def translate(genome1, genome1b, genome2, genome2b):
	trans = {}
	for gene in genome1:
		if gene.chromosome in genome1b:
			lpos = genome2.getPosition(gene.names)
			tmp = [genome2.lstGenes[pos.chromosome][pos.index].names[0] for pos in lpos if pos.chromosome in genome2b]
			if len(tmp) > 0:
				trans[gene.names[0]] = tmp
	return trans

translateMA = {}
for (x,y) in translate(genome, newGenome, ancGenome, newAncGenome).iteritems():
	assert len(y) == 1
	translateMA[x] = y[0]
print >> sys.stderr, "M->A", len(translateMA), translateMA.items()[0]

translateAM = translate(ancGenome, newAncGenome, genome, newGenome)
print >> sys.stderr, "A->M", len(translateAM), translateAM.items()[0]

# Egalite des jeux de genes consideres
assert set(translateMA.itervalues()) == set(translateAM)
assert set(translateMA) == set(itertools.chain(*translateAM.itervalues()))

def removeNewSingletons(genome, translate):
	newGenome = {}
	for chrom in genome:
		tmp = [(g,s) for (g,s) in genome[chrom] if g in translate]
		if len(tmp) >= abs(arguments["minimalLength"]):
			newGenome[chrom] = genome[chrom]
	return newGenome

while arguments["minimalLength"] < -1:
	(n1,n2) = (len(newGenome), len(newAncGenome))
	print >> sys.stderr, "iter"
	newGenome = removeNewSingletons(newGenome, translateMA)
	print >> sys.stderr, "modernGenome", sum([len(x) for x in newGenome.itervalues()])
	newAncGenome = removeNewSingletons(newAncGenome, translateAM)
	print >> sys.stderr, "ancGenome", sum([len(x) for x in newAncGenome.itervalues()])

	if (n1,n2) == (len(newGenome), len(newAncGenome)):
		print >> sys.stderr, "stop"
		break

	translateMA = {}
	for (x,y) in translate(genome, newGenome, ancGenome, newAncGenome).iteritems():
		assert len(y) == 1
		translateMA[x] = y[0]
	print >> sys.stderr, "M->A", len(translateMA), translateMA.items()[0]

	translateAM = translate(ancGenome, newAncGenome, genome, newGenome)
	print >> sys.stderr, "A->M", len(translateAM), translateAM.items()[0]

	# Egalite des jeux de genes consideres
	assert set(translateMA.itervalues()) == set(translateAM)
	assert set(translateMA) == set(itertools.chain(*translateAM.itervalues()))




# Liste des intervalles
#     - Tous
#     - Ceux avec les deux genes presents dans l'autre genome -> Eviction des genes specifiques
def listInterv(genome, translate):

	listIntAll = set()
	listIntFilt = []
	dicPos = {}
	dicLengths = {}
	for chrom in genome:
		listIntAll.update(utils.myTools.myIterator.slidingTuple(genome[chrom]))
		tmp = [(g,s) for (g,s) in genome[chrom] if g in translate]
		for (i,(g,s)) in enumerate(tmp):
			dicPos[g] = (chrom,i,s)
		dicLengths[chrom] = len(tmp)
		listIntFilt.extend(utils.myTools.myIterator.slidingTuple(tmp))
	print >> sys.stderr, len(listIntAll), len(listIntFilt), list(listIntAll)[0], listIntFilt[0]
	return (listIntAll,listIntFilt,dicPos,dicLengths)


print >> sys.stderr, "intMod",
(listIntMall,listIntMfilt,dicPosMod,dicModLengths) = listInterv(newGenome, translateMA)

print >> sys.stderr, "intAnc",
(listIntAall,listIntAfilt,dicPosAnc,dicAncLengths) = listInterv(newAncGenome, translateAM)
listIntAfilt = set(listIntAfilt)
listIntAalls = set((g1,g2) for ((g1,s1),(g2,s2)) in listIntAall)
listIntAfilts = set((g1,g2) for ((g1,s1),(g2,s2)) in listIntAfilt)
assert len(listIntAalls) == len(listIntAall)
assert len(listIntAfilts) == len(listIntAfilt)

allendsT = ["NO_END", "ONE_END", "TWO_ENDS"]

def getRoom(length, pos, after):
	return (pos == (length-1)) if after else (pos == 0)


# Parcours du point de vue de l'espece moderne
seen = set()
for ((g1,s1),(g2,s2)) in listIntMfilt:
	tg1 = translateMA[g1]
	tg2 = translateMA[g2]
	flags = ["WITHOUT_GENE_GAIN_INSIDE" if ((g1,s1),(g2,s2)) in listIntMall else "WITH_GENE_GAIN_INSIDE"]
	
	if ((tg1,s1),(tg2,s2)) in listIntAfilt:
		status = "="
		flags.append("SAME_ORIENT")
		flags.append("WITHOUT_GENE_LOSS_INSIDE" if ((tg1,s1),(tg2,s2)) in listIntAall else "WITH_GENE_LOSS_INSIDE")
		seen.add( (tg1,tg2) )
	elif ((tg2,-s2),(tg1,-s1)) in listIntAfilt:
		status = "="
		flags.append("SAME_ORIENT")
		flags.append("WITHOUT_GENE_LOSS_INSIDE" if ((tg2,-s2),(tg1,-s1)) in listIntAall else "WITH_GENE_LOSS_INSIDE")
		seen.add( (tg2,tg1) )
	
	elif (tg1,tg2) in listIntAfilts:
		status = "="
		flags.append("DIFF_ORIENT")
		flags.append("WITHOUT_GENE_LOSS_INSIDE" if (tg1,tg2) in listIntAalls else "WITH_GENE_LOSS_INSIDE")
		seen.add( (tg1,tg2) )
	elif (tg2,tg1) in listIntAfilts:
		status = "="
		flags.append("DIFF_ORIENT")
		flags.append("WITHOUT_GENE_LOSS_INSIDE" if (tg2,tg1) in listIntAalls else "WITH_GENE_LOSS_INSIDE")
		seen.add( (tg2,tg1) )

	elif tg1 == tg2:
		status = "+"
		if s1 == s2:
			flags.append("DUPLICATES_SAME_ORIENT")
		else:
			flags.append("DUPLICATES_DIFF_ORIENT")

	else:
		status = "+"
		(ac1,ai1,as1) = dicPosAnc[tg1]
		(ac2,ai2,as2) = dicPosAnc[tg2]

		room1 = getRoom(dicAncLengths[ac1], ai1, as1 == s1)
		room2 = getRoom(dicAncLengths[ac2], ai2, as2 != s2)
		flags.append(allendsT[room1 + room2])

	print "\t".join([status, "%s/%d" % (g1,s1), "%s/%d" % (g2,s2), tg1, tg2] + flags)

# Parcours du point de vue de l'ancetre (ses intervalles non conserves)
for ((g1,s1),(g2,s2)) in listIntAfilt:
	if (g1,g2) in seen:
		continue
	status = "-"
	flags = ["WITHOUT_GENE_LOSS_INSIDE" if ((g1,s1),(g2,s2)) in listIntAall else "WITH_GENE_LOSS_INSIDE"]

	allends = []
	for (tg1,tg2) in itertools.product(translateAM[g1],translateAM[g2]):
		
		(ac1,ai1,as1) = dicPosMod[tg1]
		(ac2,ai2,as2) = dicPosMod[tg2]

		room1 = getRoom(dicModLengths[ac1], ai1, as1 == s1)
		room2 = getRoom(dicModLengths[ac2], ai2, as2 != s2)
		allends.append((room1 + room2, tg1, tg2))
	if len(allends) == 0:
		print >> sys.stderr, ((g1,s1),(g2,s2)), translateAM[g1], translateAM[g2]
	(i,tg1,tg2) = max(allends)
	flags.append(allendsT[i])

	print "\t".join([status, "%s/%d" % (g1,s1), "%s/%d" % (g2,s2), tg1, tg2] + flags)


