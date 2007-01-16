#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des communautes
#

import os
import sys
import myMaths
import myGenomes
import myTools


#anc=Percomorph
#num=1
#alpha=0.48397
#relevance=0.54888
#python /users/ldog/muffato/work/scripts/extract_clusters.py $anc/log.$num $alpha $anc/nodes.$num
#| ~/work/scripts/toto2.py ~/work/phylTree.noLowCoverage.conf ~/work/data41/diags/diags.noKeepOnlyOrthos.noLowCoverage.projected.longest.moreSpecies.bz2 -ancestr=$anc
#| ~/work/scripts/printAncestralGenome.py ~/work/data41/ancGenes/ancGenes.$anc.list.bz2 /dev/stdin
#| sort -k1,1n -k2 > Genome.$anc-$alpha-$relevance

def launchCommunitiesBuild(nbItems, scoreFunc, keepLonelyNodes):
	combin = myTools.myCombinator([])
	dicAretes = dict([(i,{}) for i in xrange(nbItems)])
	for (i1,i2) in myTools.myMatrixIterator(nbItems, nbItems, myTools.myMatrixIterator.StrictUpperMatrix):
	
		score = scoreFunc(i1, i2)
	
		if score > 0:
			combin.addLink([i1,i2])
			dicAretes[i1][i2] = score
			dicAretes[i2][i1] = score

	# On traite chaque composante connexe
	res = []
	indComp = 1
	for g in combin:
	
		print >> sys.stderr, "Traitement de la composante connexe %d ..." % indComp,
		
		(stdin,stdout) = os.popen2('/users/ldog/muffato/work/scripts/walktrap/walktrap -s')
		
		nb = len(g)
		for (i1,i2) in myTools.myMatrixIterator(nb, nb, myTools.myMatrixIterator.StrictUpperMatrix):
			if g[i2] in dicAretes[g[i1]]:
				print >> stdin, i1, i2, dicAretes[g[i1]][g[i2]]
		stdin.close()

		print >> sys.stderr, "Extraction des communautes ...",

		comm = loadCommunitiesFile(stdout)
		res.extend([[g[i] for i in c] for c in comm])
		#if keepLonelyNodes:
		
		print >> sys.stderr, len(comm), "OK"


	return res

def loadCommunitiesFile(file):


	class Merge(object):
		scale = 0.0
		father = 0
		children = []

	allMerges = []

	
	def getAllChildren(father):
		alreadySeen.add(father)
		found = 0
		res = []
		for merge in allMerges:
			if merge.father == father:
				found = 1
				for child in merge.children:
					res.extend( getAllChildren(child) )
		if found == 0:
			res.append(father)
		
		return res

	
	
	for line in file:
		if line == "\n":
			break
		
		list = line.split (':')
		
		merge = Merge ()
		merge.scale = float(list[0])

		list = list[1].split('-->')
		merge.father = int(list[1])
		merge.children = [int(x) for x in list[0].split('+')]

		allMerges.append(merge)

	allMerges.reverse()

	vals = []
	for line in file:
		try:
			c = line.split()
			alpha = float(c[0])
			relevance = float(c[1])
			vals.append( (relevance, alpha) )
		except Exception:
			pass
	
	if len(vals) == 0:
		return []
	
	vals.sort()
	scale = vals[-1][1]
	print >> sys.stderr, "relevance=%f alpha=%f" % vals[-1],
	
	alreadySeen = set([])

	lstClusters = []
	for merge in allMerges:
		if merge.scale < scale and merge.father not in alreadySeen:
			sys.stderr.write('.')
			lstClusters.append( getAllChildren(merge.father) )

	if printLonelyNodes:
		for i in xrange (max / 2 + 1):
			if i not in alreadySeen:
				pass
		#		print i

	return lstClusters




