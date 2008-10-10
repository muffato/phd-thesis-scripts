#! /users/ldog/muffato/python

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""

import sys
import utils.myDiags
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

genomes = {}
ancGenes = {}

for l in sys.stdin:
	l = l.replace('\n', '')
	lstDiags = utils.myDiags.loadDiagsFile(l, phylTree.listAncestr, phylTree.officialName)

	allOK = 0.
	allPerfect = 0.
	allPerfectPairs = 0.
	allDiags = 0.
	allPairs = 0.
	allShift = 0.
	allCov = 0.

	for anc in lstDiags:
		if len(lstDiags[anc]) == 0:
			continue
		if anc not in genomes:
			genomes[anc] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[anc])
			ancGenes[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
		nbOK = 0.
		nbPerfect = 0.
		nbPerfectPairs = 0.
		nbTotPairs = 0.
		averageShift = 0.
		allPos = set()
		for (d,_,_,_,_) in lstDiags[anc]:
			lstPos = [genomes[anc].getPosition(ancGenes[anc].lstGenes[None][i].names) for i in d]
			tmp = list(set([len(x) for x in lstPos]))
			if tmp != [1]:
				print >> sys.stderr, "PB1 !!", d, lstPos, tmp
				
			lstPos = [x.pop() for x in lstPos]
			allPos.update(lstPos)
			if len(set([c for (c,_) in lstPos])) == 1:
				nbOK += 1.
				lstPos.sort()
				(_,x1) = lstPos[0]
				(_,x2) = lstPos[-1]
				if abs(x2-x1) == (len(d)-1):
					nbPerfect += 1.
			for i in xrange(len(lstPos)-1):
				(c1,x1) = lstPos[i]
				(c2,x2) = lstPos[i+1]
				if c1 != c2:
					continue
				if abs(x1-x2) == 1:
					nbPerfectPairs += 1.
				averageShift += abs(x1-x2)
				nbTotPairs += 1.

	
		# Fichier, ancetre, anciennete, %memeChr, %parfaite, ecart moyen, couverture
		print "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (l,anc, phylTree.ages[anc],100.*nbOK/len(lstDiags[anc]),100.*nbPerfect/len(lstDiags[anc]),averageShift/nbTotPairs,100.*nbPerfectPairs/nbTotPairs,100. * float(len(allPos)) / float(sum([len(x) for x in genomes[anc].lstGenes.itervalues()])))
		
		allOK += nbOK
		allPerfect += nbPerfect
		allPerfectPairs += nbPerfectPairs
		allDiags += len(lstDiags[anc])
		allPairs += nbTotPairs
		allShift += averageShift
		allCov += float(len(allPos)) / float(sum([len(x) for x in genomes[anc].lstGenes.itervalues()]))


	print  "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (l,"----",100.*allOK/allDiags, 100.*allPerfect/allDiags, allShift/allPairs, 100.*allPerfectPairs/allPairs, 100.*allCov/len(lstDiags))

