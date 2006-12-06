#! /users/ldog/muffato/python

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags


#############
# FONCTIONS #
#############

def loadDiagsFile(nom, diagEntry, genesAnc):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	for l in f:
		ct = l.split('\t')
		anc = ct[0]
		d1 = ct[4].split()
		d2 = ct[7].split()
		d = d1+d2
		da1 = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d1]
		da2 = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d2]
		if "" in da1:
			da = da2
		else:
			da = da1
		diagEntry[anc].append( (d, da, (ct[2],ct[3],d1,da1), (ct[5],ct[6],d2,da2), int(ct[1])) )
		
	f.close()
	print >> sys.stderr, "OK" #, " ".join(["%s:%d" % (anc,len(diagEntry[anc])) for anc in diagEntry])
	


def extractLongestOverlappingDiags(oldDiags, newDiags):

	dic = {}
	for i in xrange(len(oldDiags)):
		d = oldDiags[i][0]
		for s in d:
			if s not in dic:
				dic[s] = []
			dic[s].append(i)
	
	combin = utils.myTools.myCombinator([[x] for x in xrange(len(oldDiags))])
	for s in dic:
		combin.addLink(dic[s])

	for g in combin:
		for res in utils.myDiags.getLongestPath([oldDiags[i][1] for i in g]):
			ok = set([])
			for i in g:
				d = oldDiags[i][1]
				flag = False
				for j in xrange(len(d)-1):
					if (d[j] not in res[0]) or (d[j+1] not in res[0]):
						continue
					if abs(res[0].index(d[j])-res[0].index(d[j+1])) == 1:
						flag = True
						break
				if flag:
					ok.add( (oldDiags[i][2][0],oldDiags[i][2][1]) )
					ok.add( (oldDiags[i][3][0],oldDiags[i][3][1]) )
			print "%s\t%d\t%s\t%s" % (anc, len(res[0]), " ".join([str(x) for x in res[0]]), " ".join(["%s.%s" % (e,c) for (e,c) in ok]))
			newDiags.append( (len(res[0]), res[0], ok) )
	

def combinSameChr():
	combin = utils.myTools.myCombinator([])
	#fils = [['Human'], ['Chimp']]
	fils = phylTree.getBranchesSpecies(sys.argv[3])
	for (i,j) in utils.myTools.myMatrixIterator(len(diagFinal), len(diagFinal), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		commun = set([e for (e,_) in diagFinal[i][2].intersection(diagFinal[j][2])])
		diff = set([e for (e,_) in diagFinal[i][2].symmetric_difference(diagFinal[j][2])])
		filsOK = [len(commun.intersection(x)) for x in fils]
		filsNO = [len(diff.intersection(x)) for x in fils]
		outgroupOK = len(commun) - sum(filsOK)
		outgroupNO = len(diff) - sum(filsNO)
		#print diagFinal[i][2], diagFinal[j][2]
		#print commun, diff, filsOK, filsNO, outgroupOK, outgroupNO
		#if max(filsNO) == 0:
		#if (min(filsOK) >= 1) or (max(filsOK) >= 1 and outgroupOK >= 1):
		if min(filsOK) >= 1:
		#if outgroupOK >= 1 and min(filsOK) >= 1 and max(filsOK) >= 1:
			combin.addLink([i,j])
	for g in combin:
		r = set(utils.myMaths.flatten([diagFinal[i][1] for i in g]))
		print " ".join([str(x) for x in r])


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf", "diagsList"], \
	[("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
#geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)

# Pour sauver de la memoire
#del geneBank.dicEspeces
#for esp in geneBank.dicEspeces:
#	del geneBank.dicEspeces[esp].dicGenes

diagEntry = {}
genesAnc = {}
for anc in phylTree.items:
	diagEntry[anc] = []
	genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc, False, False)

loadDiagsFile(noms_fichiers["diagsList"], diagEntry, genesAnc)

# 2. On extrait les diagonales les plus longues
for anc in phylTree.items:
	print >> sys.stderr, "Extraction des chevauchements les plus longs de %s ..." % anc,
	tmp = []
	extractLongestOverlappingDiags(diagEntry[anc], tmp)
	print >> sys.stderr, "OK (%d -> %d)" % (len(diagEntry[anc]), len(tmp))
	diagEntry[anc] = tmp
	

