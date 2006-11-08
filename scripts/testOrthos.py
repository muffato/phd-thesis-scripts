#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales de genes entre deux especes
"""

##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags
import cPickle

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesAnc"], \
	[("ancetre",str,""), ("output",str,""), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

def f1():
	tabGenesAnc = []

	for nom in sys.argv[2:]:
		tmp = utils.myGenomes.loadGenome(nom)
		del tmp.lstGenes
		tabGenesAnc.append(tmp)

	genomeAnc = utils.myGenomes.loadGenome(sys.argv[1])
	del genomeAnc.dicGenes
	print "\t\t%s" % "\t".join([str(x) for x in genomeAnc.lstChr])
	for g in genomeAnc:
		score = dict([(c,0) for c in genomeAnc.lstChr])
		for gen in tabGenesAnc:
			s = g.names[0]
			if s in gen.dicGenes:
				score[gen.dicGenes[s][0]] += 1
		print "%s\t%d\t%s" % (g.chromosome, g.beginning, "\t".join([str(score[x]) for x in score]))



def f2():
	f = open(noms_fichiers["genesAnc"], 'r')

	lst = f.readlines()[7:-1]
	f.close()

	for i in range(len(lst)):
		continue
		c = lst[i].split()
		for j in range(len(c)):
			if float(c[j]) > 0.5 and float(c[j]) < 1.5:
				print i,j+i+1
			



genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genesAnc"])

nb = 1
for l in sys.stdin:
	c = l.split()
	for i in c:
		print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][int(i)].names)
	nb += 1

sys.exit(0)
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
nbEsp = len(phylTree.getSpecies(phylTree.root))

for anc in phylTree.items:
	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	fils = utils.myMaths.flatten(groupes)

	nbO = nbEsp-len(fils)
	nbA = len(groupes[0])
	nbB = len(groupes[1])
	print anc, nbA*nbB + nbA*nbO + nbB*nbO
