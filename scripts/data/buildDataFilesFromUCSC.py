#! /users/ldog/muffato/python -OO


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import operator
import utils.myPhylTree
import utils.myTools


#############
# FONCTIONS #
#############


def FUNCnet():

	links = {}

	lizChr = None
	f = utils.myTools.myOpenFile(noms_fichiers["chains"], "r")

	for l in f:
		c = l.split()
		if "net" in l:
			lizChr = c[1]
			continue
		if 'fill' not in l:
			continue
		#lizChr = c[3]
		liz1 = int(c[1])
		lizSize = int(c[2])

		lizGt = list(genomeLizard.getGenesAt(lizChr, liz1, liz1+lizSize))
		if len(lizGt) > 0:
			print lizGt
			print
	f.close()


def FUNCnetaxt():

	nb = set()
	f = utils.myTools.myOpenFile(noms_fichiers["chains"], "r")
	for l in f:
		c = l.split()
		lizChr = c[1]
		liz1 = int(c[2])
		liz2 = int(c[3])
		chickChr = utils.myGenomes.commonChrName(c[4][3:])
		chickStrand = c[7]
		chick1 = int(c[5])
		chick2 = int(c[6])
		#lizGt = list(genomeLizard.getGenesAt(lizChr, liz1, liz2))
		#if len(lizGt) > 0:
		#	print lizGt
		#	print
		if chickStrand == "+":
			chickGt = list(genomeChicken.getGenesAt(chickChr, chick1, chick2))
			nb.update(chickGt)
			#if len(chickGt) > 0:
			#	print chickGt
			#	print
	f.close()
	print >> sys.stderr, len(nb)



def FUNCallchain(genome, chain):

	genomeChicken = utils.myGenomes.EnsemblGenome(genome)
	f = utils.myTools.myOpenFile(chain, "r")
	for l in f:
		if not l.startswith("chain"):
			continue
		c = l.split()
		lizChr = c[2]
		lizStrand = c[4]
		liz1 = int(c[5])
		liz2 = int(c[6])
		chickChr = utils.myGenomes.commonChrName(c[7][3:])
		chickStrand = c[9]
		chick1 = int(c[10])
		chick2 = int(c[11])
		if lizStrand == "-":
			# On inverse les coordonnees du lezard
			lizChrSize = int(c[3])
			(liz1,liz2) = (lizChrSize-1-liz2,lizChrSize-1-liz1)
		if chickStrand == "-":
			# On inverse les coordonnees poulet
			chickChrSize = int(c[8])
			(chick1,chick2) = (chickChrSize-1-chick2,chickChrSize-1-chick1)
		chickGt = list(genomeChicken.getGenesAt(chickChr, chick1, chick2))
		lizGt = list(genomeLizard.getGenesAt(lizChr, liz1, liz2))
		tmp = [g.names[0] for g in lizGt+chickGt]
		if len(chickGt) > 0 and len(lizGt) > 0:
			print " ".join(tmp)
		#for g1 in chickGt:
		#	for g2 in lizGt:
		#		print
		#		linksCL.setdefault(g1,set()).add(g2)
		#		linksLC.setdefault(g2,set()).add(g1)
			
	f.close()
	#for g1 in linksCL:
	#	if len(linksCL[g1]) == 1:
	#		g2 = linksCL[g1].pop()
	#		if len(linksLC[g2]) == 1:
	#			print g1.names[0], g2.names[0]
	#	print g1.names[0], " ".join([g.names[0] for g in linksCL[g1]])
	#for g1 in linksLC:
	#	print g1.names[0], " ".join([g.names[0] for g in linksLC[g1]])



# Arguments
genomeLizard = utils.myGenomes.EnsemblGenome(sys.argv[1])

nbEntitites = len(sys.argv)/2 - 1
for i in xrange(nbEntitites):
	FUNCallchain(sys.argv[2*i+2], sys.argv[2*i+3])


#FUNCnetaxt()
#FUNCallchain()


