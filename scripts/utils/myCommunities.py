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

	print >> sys.stderr, "Calcul des scores des aretes ...",
	combin = myTools.myCombinator([])
	dicAretes = dict([(i,{}) for i in xrange(nbItems)])
	for (i1,i2) in myTools.myMatrixIterator(nbItems, nbItems, myTools.myMatrixIterator.StrictUpperMatrix):
	
		score = scoreFunc(i1, i2)
		if score > 0:
			combin.addLink([i1,i2])
			dicAretes[i1][i2] = score

	print >> sys.stderr, "Division en communautes ..."
	res = []
	for g in combin:
		tmp = launchCommunitiesBuild1(g, dicAretes)
		if len(tmp) == 0:
			tmp = [(1,0,[g],[])]
		res.append(tmp)
	print >> sys.stderr, "OK"
	
	return res



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
		#print >> sys.stderr, "New result alpha=%f relevance=%f %d clusters %d/%d N/A" % (alpha, relevance, len(newC), len(lonelyNodes), nb)
		
	return res



def loadCommunitiesFileB(file):

	class Merge(object):
		scale = 0.0
		father = 0
		children = []
	
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

	
	allMerges = []
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

