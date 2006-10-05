#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["geneList.conf", "genomeMammals", "genomeOutgroup", "orthologuesList"], \
	[], \
	"Extrait tous les groupes de genes qui sont sur un meme chromosome de maniere tres stringente" \
)


# 1. On lit tous les fichiers
geneBank = myOrthos.GeneBank(noms_fichiers[0], "OC")
genomeMammals = myOrthos.AncestralGenome(noms_fichiers[1], True)
genomeOutgroup = myOrthos.AncestralGenome(noms_fichiers[2], True)
genesAnc = myOrthos.AncestralGenome(noms_fichiers[3], False)

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
	for (i,j) in myTools.myMatrixIterator(len(lstGenes), len(lstGenes), myTools.myMatrixIterator.StrictUpperMatrix):
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

lstGenes = genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr]
geneBank.dicEspeces['T'] = genomeOutgroup
geneBank.dicEspeces['A'] = genomeMammals
buildAncestrGenome(geneBank.dicEspeces, lstGenes)

