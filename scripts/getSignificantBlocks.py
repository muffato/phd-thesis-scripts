#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
Trie les gens selon l'ordre consensus issu des especes de l'arbre phylogenetique
"""

import sys
import math
import itertools
import collections

import utils.myGenomes
import utils.myTools
import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("ancGenesFile",str,""), ("diagsFile",str,""), ("nbSamples",int,1000), ("pvalue",float,2.5), ("withSingletons",bool,False)], __doc__ )


phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

limsign = int(arguments["nbSamples"] * arguments["pvalue"] / 100.)


def do(anc):
	# Lecture des tailles des blocs de syntenie
	f = utils.myFile.myTSV.reader(arguments["diagsFile"] % phylTree.fileName[anc])
	tailles = []
	for l in f[1]:
		tailles.append(int(l[1]))
	f[0].close()

	# Eventuellement les singletons
	if arguments["withSingletons"]:
		singletons = len(utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc]).lstGenes[None]) - sum(tailles)
		tailles.extend([1] * singletons)
	
	print anc

	count = utils.myMaths.myStats.count(tailles)
	#print count

	# Les parametres du modele geometrique
	m = utils.myMaths.myStats.mean(tailles)
	lp = [1/m]
	stddev = utils.myMaths.myStats.stddev(tailles, m)
	if stddev > 0:
		variance = stddev**2
		lp.append((math.sqrt(1+4*variance)-1)/(2*variance))
		lp.append(sum(lp)/2)
	
	# Test de chaque p
	for p in lp:
		
		print "**** with p=%.3f ****" % p
		
		# La queue de la distribution a 10^-X
		for alpha in itertools.count(2):
			lim = int(1-(alpha*math.log(10))/math.log(1-p))
			n = len([l for l in tailles if l >=lim])
			print "%d blocks with pvalue < 10^-%d (min-length=%d)" % (n, alpha, lim)
			if n == 0:
				break
		print

		# C'est parti pour les tirages aleatoires
		allcounts = []
		m = max(count)
		# Le bon nombre d'echantillons
		for _ in xrange(arguments["nbSamples"]):
			# Autant de tailles qu'observe
			newcount = collections.defaultdict(int)
			for _ in xrange(len(tailles)):
				newcount[utils.myMaths.randomValue.geometric(p)] += 1
			allcounts.append(newcount)
			m = max(m, max(newcount))

		for l in xrange(1, m+1):
			seen = [c[l] for c in allcounts]
			seen.sort()
			# Moins qu'attendu
			if (count[l] < seen[limsign]):
				print "diff on size %d: observed=%d << min-expected=%d" % (l,count[l],seen[limsign])
			# Plus qu'attendu
			if (count[l] > seen[-limsign-1]):
				print "diff on size %d: observed=%d >> max-expected=%d" % (l,count[l],seen[-limsign-1])
		print


for anc in phylTree.listAncestr:
	do(anc)

