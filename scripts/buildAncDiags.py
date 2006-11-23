#! /users/ldog/muffato/python

__doc__ = """
Extrait toutes les diagonales entre chaque paire d'especes.
Les diagonales apportent les genes qui etaient sur un meme chromosome
  depuis leur ancetre commun dans les deux lignees.
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

def loadAncGenesFile(anc, genomes, locations):

	# Le fichier de genes ancestraux
	genesAnc = utils.myGenomes.loadGenome(options["ancGenesFile"] % anc)
	nbGenesAnc = len(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr])
	del genesAnc.lstGenes

	# Les listes des especes entre lesquelles on cherche des diagonales
	groupes = phylTree.getBranchesSpecies(anc)
	fils = utils.myMaths.flatten(groupes)
	
	# Traduction des genomes en liste des genes ancestraux
	print >> sys.stderr, "Distribution des genes de %s " % anc,
	for e in fils:
		newLoc = [[] for x in xrange(nbGenesAnc)]
		genome = geneBank.dicEspeces[e]
		newGenome = {}
		for c in genome.lstChr:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
			if not options["keepOrthosLess"]:
				newGenome[c] = [x for x in newGenome[c] if x[0] != -1]
			for i in xrange(len(newGenome[c])):
				(ianc,s) = newGenome[c][i]
				if ianc != -1:
					newLoc[ianc].append( (c,i,s) )
		genomes[e] = newGenome
		locations[e] = newLoc
		sys.stderr.write(".")
	
	print >> sys.stderr, " OK"


def calcDiags():

	# La fonction qui permet de traiter les diagonales
	def combinDiag(c1, c2, d1, d2):
		global diagEntry

		if len(d1) < options["minimalLength"]:
			return
		
		for (e, tmp) in toStudy:
			if e == e1:
				dd = [genomes[tmp][e][c1][i][0] for i in d1]
			else:
				dd = [genomes[tmp][e][c2][i][0] for i in d2]
			
			diagEntry[tmp].addDiag(dd, [(e1,c1), (e2,c2)] )
			
	
	n = max([len(x) for x in listEspeces])
	for i in xrange(len(listEspeces)):
		print >> sys.stderr, ("Extraction des diagonales de %-"+str(n) + "s ") % listEspeces[i],
		for j in xrange(len(listEspeces)):
			if j < i:
				print >> sys.stderr, '-',
				continue
			elif j == i:
				print >> sys.stderr, 'X',
				continue
				
			e1 = listEspeces[i]
			e2 = listEspeces[j]
			anc = phylTree.getFirstParent(e1,e2)
			toStudy = [(e1,anc)]
			for tmp in phylTree.items:
				s = phylTree.getSpecies(tmp)
				if (e1 in s) and (e2 not in s):
					toStudy.append( (e1,tmp) )
				elif (e2 in s) and (e1 not in s):
					toStudy.append( (e2,tmp) )
				
			utils.myDiags.iterateDiags(genomes[anc][e1], locations[anc][e2], options["fusionThreshold"], options["sameStrand"], combinDiag)
			print >> sys.stderr, '+',
		
		for anc in phylTree.items:
			if e1 in genomes[anc]:
				del genomes[anc][e1]
				del locations[anc][e1]
		
		print >> sys.stderr, "OK"

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOrthosLess",bool,False), \
	("checkInsertions",bool,False), ("extendLeftRight",bool,False), ("minOverlap",float,-1), ("combinSameChr",bool,False), ("checkCliques",int,-1), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)

# Pour sauver de la memoire
for esp in geneBank.dicEspeces:
	del geneBank.dicEspeces[esp].dicGenes

# 2. On prepare tous les genomes ancestraux, les genomes traduits ...
genomes = dict( [(anc, {}) for anc in phylTree.items] )
locations = dict( [(anc, {}) for anc in phylTree.items] )
for anc in phylTree.items:
	loadAncGenesFile(anc, genomes[anc], locations[anc])

# Plus besoin de ca ...
del geneBank

# La structure qui accueillera les diagonales et le calcul
diagEntry = dict( [(anc, utils.myDiags.DiagRepository()) for anc in phylTree.items] )
calcDiags()

# Plus besoin de ca ...
del genomes
del locations

for anc in diagEntry:
	
	print >> sys.stderr, "Traitement de %s ..." % anc,
	lst = diagEntry[anc]
	print >> sys.stderr, lst.nbRealDiags(),
	
	if options["checkInsertions"]:
		it = 0
		while lst.checkInsert():
			it += 1
		print >> sys.stderr, "I [%d]" % it, lst.nbRealDiags(),
	
	if options["extendLeftRight"]:
		lst.buildVoisins()
		it = 0
		while (lst.extendLeft() or lst.extendRight()):
			it += 1
		print >> sys.stderr, "E [%d]" % it, lst.nbRealDiags(),
	
	if options["checkCliques"] > 0:
		
		nb = 0
		it = 0
		while nb != lst.nbRealDiags():
			nb = lst.nbRealDiags()
			lst.buildCliques()
			if len(lst.cliquesList) <= options["checkCliques"]:
				continue
			cl = lst.cliquesList[-1]
			if len(cl) == 0:
				continue
			s = set([])
			for i in cl.pop():
				s.update(lst.lstDiagsSet[i])
			lst.addDiag(list(s), [])
			it += 1

		print >> sys.stderr, "L [%d]" % it, lst.nbRealDiags(),

	if options["minOverlap"] > 0:
		lst = lst.combinOverlap(options["minOverlap"])
		print >> sys.stderr, "O", lst.nbRealDiags(),

	if options["combinSameChr"]:
		lst.combinDiags(phylTree.getBranchesSpecies(anc))
		print >> sys.stderr, "C", lst.nbRealDiags(),


	for (d,_,l) in lst:
		#print anc, " ".join([str(x) for x in d])
		print "%s\t%s\t%s" % (anc, " ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l]))
	
	print >> sys.stderr, "OK"

