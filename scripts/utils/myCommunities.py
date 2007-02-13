#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des communautes
#

import os
import sys
import myMaths
import myGenomes
import myTools


def launchCommunitiesBuildB(nbItems, scoreFunc):
	
	combin = myTools.myCombinator([])
	dicAretes = dict([(i,{}) for i in xrange(nbItems)])
	for (i1,i2) in myTools.myMatrixIterator(nbItems, nbItems, myTools.myMatrixIterator.StrictUpperMatrix):
	
		score = scoreFunc(i1, i2)
		if score > 0:
			combin.addLink([i1,i2])
			dicAretes[i1][i2] = score

	res = []
	for g in combin:
		#tmp = launchCommunitiesBuild1(g, dicAretes)
		#if len(tmp) == 0:
		tmp = [(1,0,[g],[])]
		res.append(tmp)
	return res
	#return [launchCommunitiesBuild1(g, dicAretes) for g in combin]



def launchCommunitiesBuild1(nodes, edges):
	
	indNodes = {}
	nb = len(nodes)
	for i in xrange(nb):
		indNodes[nodes[i]] = i
	
	(stdin,stdout) = os.popen2('/users/ldog/muffato/work/scripts/walktrap -s -t10')
	
	for x in nodes:
		if x not in edges:
			continue
		for y in nodes:
			if y in edges[x]:
				print >> stdin, indNodes[x], indNodes[y], edges[x][y]
	stdin.close()

	res = []
	for (alpha,relevance,clusters) in loadCommunitiesFileB(stdout):
		newC = [[nodes[i] for i in c] for c in clusters]
		lonelyNodes = set(nodes).difference(myMaths.flatten(newC))
		res.append( (alpha, relevance, newC, list(lonelyNodes)) )
		print >> sys.stderr, "New result alpha=%f relevance=%f %d clusters %.2f N/A" % \
		(alpha, relevance, len(newC), (100*len(lonelyNodes))/nb)
		
	return res



def loadCommunitiesFileB(file):


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

	res = []
	for line in file:
		print >> sys.stderr, line,
		try:
			c = line.split()
			scale = float(c[0])
			relevance = float(c[1])

			alreadySeen = set([])
			lstClusters = []
			for merge in allMerges:
				if merge.scale < scale and merge.father not in alreadySeen:
					lstClusters.append( getAllChildren(merge.father) )

			res.append( (scale, relevance, lstClusters) )

		except Exception:
			pass

	return res

















def launchCommunitiesBuild(nbItems, scoreFunc, keepLonelyNodes = False, minRelevance = 0, minCoverage = 0, bestRelevance = True):
	
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
		
		(stdin,stdout) = os.popen2('/users/ldog/muffato/work/scripts/walktrap -s -t10')
		
		nb = len(g)
		for (i1,i2) in myTools.myMatrixIterator(nb, nb, myTools.myMatrixIterator.StrictUpperMatrix):
			if g[i2] in dicAretes[g[i1]]:
				print >> stdin, i1, i2, dicAretes[g[i1]][g[i2]]
		stdin.close()

		print >> sys.stderr, "Extraction des communautes ...",

		(v, comm) = loadCommunitiesFile(stdout, bestRelevance)
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

def loadCommunitiesFile(file, bestRelevance):


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
		print >> sys.stderr, line,
		try:
			c = line.split()
			scale = float(c[0])
			relevance = float(c[1])
			if bestRelevance:
				vals.append( (relevance, scale) )
			else:
				vals.append( (scale, relevance) )
		except Exception:
			pass
	
	if len(vals) == 0:
		return (0,[])

	vals.sort()
	if bestRelevance:
		(relevance, scale) = vals[-1]
	else:
		(scale, relevance) = vals[-1]
		
	print >> sys.stderr, "relevance=%f alpha=%f " % (relevance, scale),
	
	alreadySeen = set([])

	lstClusters = []
	for merge in allMerges:
		if merge.scale < scale and merge.father not in alreadySeen:
			sys.stderr.write('.')
			lstClusters.append( getAllChildren(merge.father) )

	return (relevance, lstClusters)




