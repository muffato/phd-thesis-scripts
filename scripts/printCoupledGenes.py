#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["geneList.conf", "orthologuesList"], \
	[], \
	"Extrait tous les groupes de genes qui sont sur un meme chromosome de maniere tres stringente" \
)


# 1. On lit tous les fichiers
genomesList = "HDMWOC"
geneBank = utils.myGenomes.GeneBank(noms_fichiers[0], genomesList)
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers[1], False)

def buildAncestrGenome(genomesList, lstGenes):

	def speciesString(s):
		return '"' + "".join([e for e in s]) + '"'

	print >> sys.stderr, "1.",
	#tabE = dict([(e,{}) for e in genomesList])
	tabAnc = [set([]) for x in lstGenes]
	for i in xrange(len(lstGenes)):
		g = lstGenes[i]
		for s in g.names:
			for e in genomesList:
				if s in genomesList[e].dicGenes:
					(c,_) = genomesList[e].dicGenes[s]
					#tabE[e][i] = c
					tabAnc[i].add( (e,c) )
					break

	print >> sys.stderr, "2.",
	for (i,j) in utils.myTools.myMatrixIterator(len(lstGenes), len(lstGenes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		#nbOK = set([])
		#nbNO = set([])
		#for e in genomesList:
		#	if (i not in tabE[e]) or (j not in tabE[e]):
		#		continue
		#	if tabE[e][i] == tabE[e][j]:
		#		nbOK.add(e)
		#	else:
		#		nbNO.add(e)
		
		nbOK = [e for (e,_) in tabAnc[i] & tabAnc[j]]
		
		nbNO = (set([e for (e,_) in tabAnc[i]]) & set([e for (e,_) in tabAnc[j]])).difference(nbOK)
		
		#if len(a) != len(nbOK) or len(b) != len(nbNO):
		#	print >> sys.stderr, tabAnc[i], tabAnc[j], a, b, nbOK, nbNO

		if len(nbOK) != 0:
			print "%d\t%d\t%s\t%s" % (i,j,speciesString(nbOK),speciesString(nbNO))

	print >> sys.stderr, "OK"

lstGenes = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
buildAncestrGenome(geneBank.dicEspeces, lstGenes)

