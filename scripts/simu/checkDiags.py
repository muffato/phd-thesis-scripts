#! /users/ldog/muffato/python -OO

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPhylTree

#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom):
	
	print >> sys.stderr, "Chargement de %s ..." % nom,
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = {}
	for l in f:
		
		# On n'utilise que "l'appartenance" a un chromosome, le fait que ce soit du random n'est pas important
		ct = l.replace('\n', '').replace("_random", "").split('\t')
		d = [int(x) for x in ct[2].split(' ')]
		esp = set()
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp if (e in phylTree.officialName) and ('Un' not in c)])
		if ct[0] in lst:
			lst[ct[0]].append( (d,esp) )
		else:
			lst[ct[0]] = [(d,esp)]

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % sum([len(lst[x]) for x in lst])
	return lst


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

genomes = {}
ancGenes = {}
for anc in phylTree.listAncestr:
	genomes[anc] = utils.myGenomes.Genome(options["genesFile"] % phylTree.fileName[anc])
	ancGenes[anc] = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])



for l in sys.stdin:
	l = l.replace('\n', '')
	lstDiags = loadDiagsFile(l)

	allOK = 0.
	allPerfect = 0.
	allPerfectPairs = 0.
	allDiags = 0.
	allPairs = 0.
	allShift = 0.
	allCov = 0.

	for anc in lstDiags:
		if anc not in genomes:
			continue
		nbOK = 0.
		nbPerfect = 0.
		nbPerfectPairs = 0.
		nbTotPairs = 0.
		averageShift = 0.
		allPos = set()
		for (d,_) in lstDiags[anc]:
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


		print "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (l,anc, phylTree.ages[anc],100.*nbOK/len(lstDiags[anc]),100.*nbPerfect/len(lstDiags[anc]),averageShift/nbTotPairs,100.*nbPerfectPairs/nbTotPairs,100. * float(len(allPos)) / float(sum([len(x) for x in genomes[anc].lstGenes.itervalues()])))
		
		allOK += nbOK
		allPerfect += nbPerfect
		allPerfectPairs += nbPerfectPairs
		allDiags += len(lstDiags[anc])
		allPairs += nbTotPairs
		allShift += averageShift
		allCov += float(len(allPos)) / float(sum([len(x) for x in genomes[anc].lstGenes.itervalues()]))


	print  "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (l,"----",100.*allOK/allDiags, 100.*allPerfect/allDiags, allShift/allPairs, 100.*allPerfectPairs/allPairs, 100.*allCov/len(lstDiags))

