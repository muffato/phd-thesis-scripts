#! /users/ldog/muffato/python

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
-> Variante pour ressembler a la vraie evolution:
    - Contraintes de rearrangements sur poulet/opossum/rongeurs
    - Especes en cours d'assemblage: orny/xenope
    - Especes a 2X de couverture
"""

import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs(
	[("phylTree.conf",file), ("root",str), ("iniChrNumber",int)], [
	
	# Taux de rearrangements
	("rate:eventMaxAccel",float,1.5), ("rate:vonMisesKappa",float,2.),

	# Parametres des genes
	("gene:statsFile",str,""), ("gene:treesSignatures",str,""), ("gene:clusteredGenesRatio",float,.2),
	("gene:initialNumber",int,22000), ("gene:dupNumber",float,1.75), ("gene:lossRate",float,19.), ("gene:gainRate",float,5.), ("gene:dupRate",float,8.),
	("gene:clusterDupLength",float,1.14), ("gene:clusterGainLength",float,1.36), ("gene:clusterLossLength",float,2.46),
	
	# Rearrangements de chromosomes
	("chr:ratesPhylTree",str,""), ("chr:rateMultiplier",float,1.),
	("chr:invertRate",float,1.), ("chr:translocRate",float,0.3), ("chr:fusionRate",float,0.05), ("chr:breakRate",float,0.05),
	("chr:vonMisesMean",float,.2), ("chr:vonMisesKappa",float,2.),
	
	# Fichiers
	("out:genomeFile",str,"simu/genes/genes.%s.list.bz2"),
	("out:ancGenesFile",str,"simu/ancGenes/ancGenes.%s.list.bz2")],
	__doc__
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

print >> sys.stderr, "Building %s gene list" % arguments["root"],
def do():
	geneStatsFile = utils.myFile.openFile(arguments["gene:statsFile"], "r")
	geneNames = eval(geneStatsFile.readline())
	
	# Rassemble les signatures identiques -> sauve de la memoire
	bysign = collections.defaultdict(list)
	f = utils.myFile.openFile(arguments["gene:treesSignatures"], "r")
	for l in f:
		(name,_,val) = l.partition(' ')
		bysign[tuple(eval(val.replace('None', '0'))[1:])].append(name)
	f.close()
	
	# On revient au tableau associatif de base
	byname = {}
	for (s,l) in bysign.iteritems():
		for x in l:
			byname[x] = s
	assert set(byname) == set(geneNames)
	
	#@utils.myTools.memoize
	def prod(v1, v2):
		s = 0.
		n = 0
		for (x,y) in itertools.izip(v1, v2):
			s += x & y
			n += x | y
		return s/n if n != 0 else -1
	
	for (s1,s2) in utils.myTools.myIterator.tupleOnUpperList(bysign.keys()):
		s = prod(s1, s2)
		if s1 is s2:
			for (g1,g2) in utils.myTools.myIterator.tupleOnUpperList(bysign[s1]):
				print g1, g2, s
			#print (len(bysign[s1])*(len(bysign[s1])-1))/2
		else:
			for (g1,g2) in itertools.product(bysign[s1], bysign[s2]):
				print g1, g2, s
				#print len(bysign[s1])*len(bysign[s2])

do()
print >> sys.stderr, "OK"

