#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des communautes
#

import os
import sys
import myMaths
import myGenomes
import myTools

def launchCommunitiesBuild(nbItems, scoreFunc, keepLonelyNodes = False, minRelevance = 0, minCoverage = 0):
	
	print >> sys.stderr, "Calcul des scores des aretes ...",
	combin = myTools.myCombinator([])
	dicAretes = dict([(i,{}) for i in xrange(nbItems)])
	for (i1,i2) in myTools.myMatrixIterator(nbItems, nbItems, myTools.myMatrixIterator.StrictUpperMatrix):
	
		score = scoreFunc(i1, i2)
	
		if score > 0:
			combin.addLink([i1,i2])
			dicAretes[i1][i2] = score
			dicAretes[i2][i1] = score
	print >> sys.stderr, "OK"

	# On traite chaque composante connexe
	res = []
	relev = []
	
	indComp = 0
	for g in combin:
	
		indComp += 1
		print >> sys.stderr, "Traitement de la composante connexe %d ..." % indComp,
		
		(stdin,stdout) = os.popen2('/users/ldog/muffato/work/scripts/walktrap/walktrap -s -t10')
		
		nb = len(g)
		for (i1,i2) in myTools.myMatrixIterator(nb, nb, myTools.myMatrixIterator.StrictUpperMatrix):
			if g[i2] in dicAretes[g[i1]]:
				print >> stdin, i1, i2, dicAretes[g[i1]][g[i2]]
		stdin.close()

		print >> sys.stderr, "Extraction des communautes ...",

		(v, comm) = loadCommunitiesFile(stdout)
		tout = set(xrange(nb))
		r = []
		for c in comm:
			tout.difference_update(c)
			r.append([g[i] for i in c])
		
		
		# Lorsque trop de diagonales sont laissees de cote, c'est qu'il
		# y a un probleme, on ne divise pas alors la composante connexe
		pourcentage = (100*len(tout))/nb
		print >> sys.stderr, " %d%% N/A," % pourcentage,

		if (pourcentage > 100*(1-minCoverage)) or (len(r) == 1) or (v < minRelevance):
			r = [g[:]]
		else:
			if keepLonelyNodes:
				r.append([g[i] for i in tout])
		
		res.extend(r)
		relev.append(v)
		print >> sys.stderr, "%d clusters" % len(r)


	return (relev, res)

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
		return (0,[])
	
	vals.sort()
	scale = vals[-1][1]
	print >> sys.stderr, "relevance=%f alpha=%f " % vals[-1],
	
	alreadySeen = set([])

	lstClusters = []
	for merge in allMerges:
		if merge.scale < scale and merge.father not in alreadySeen:
			sys.stderr.write('.')
			lstClusters.append( getAllChildren(merge.father) )

	return (vals[-1][0], lstClusters)




