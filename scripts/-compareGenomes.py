#! /users/ldog/muffato/python -OO

# Librairies
import sys
import utils.myGenomes
import utils.myTools

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["orthologuesList"], [], "Compare les genomes indiques sur l'entree standard et renvoie les stats sur la separation des paires de genes")

# Chargement des fichiers
genomes = []
for l in sys.stdin:
	for s in l.split():
		genomes.append(utils.myGenomes.Genome(s))
genesAnc = utils.myGenomes.Genome(noms_fichiers["orthologuesList"])

# Reecrit le genome avec les numeros des genes ancestraux
def rewriteGenome(genome):
	newGen = {}
	genesUsed = set()
	for (c,chrom) in genome.lstGenes.iteritems():
		pos = [genesAnc.getPosition(g.names) for g in chrom]
		newChrom = [x.pop()[1] for x in pos if len(x) >= 1]
		newGen[c] = newChrom
		genesUsed.update(newChrom)
	return (newGen,genesUsed)

genomes = [rewriteGenome(g) for g in genomes]
del genesAnc

# Recupere les genes communs aux deux genomes
genesUsed = genomes[0][1]
for i in xrange(1,len(genomes)):
	genesUsed.intersection_update(genomes[i][1])

# Fait la correspondance gene -> chromosome pour les genes en commun
def mkDic(genome):
	dic = {}
	for (c,chrom) in genome.iteritems():
		for g in genesUsed.intersection(chrom):
			dic[g] = c
	return dic

dic = [mkDic(g) for (g,_) in genomes]
del genomes
genesUsed = list(genesUsed)

nb = [0] * (2 ** len(dic))

for (g1,g2) in utils.myTools.myIterator.tupleOnStrictUpperList(genesUsed):
	x = 0
	pos = 1
	for d in dic:
		x += (d[g1] == d[g2]) * pos
		pos *= 2
	nb[x] += 1

print sum(nb)
print nb
print [(100.*i) / sum(nb) for i in nb]

