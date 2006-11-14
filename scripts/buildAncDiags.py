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
	genesAnc = utils.myGenomes.loadGenome(options["ancGenesFile"] % anc)

	# Les listes des especes entre lesquelles on cherche des diagonales
	groupes = phylTree.getBranchesSpecies(anc)
	fils = utils.myMaths.flatten(groupes)
	
	# Traduction des genomes en liste des genes ancestraux
	print >> sys.stderr, "Traduction avec les genes de", anc, "",
	for e in fils:
		genome = geneBank.dicEspeces[e]
		newGenome = {}
		for c in genome.lstChr:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
		genomes[e] = newGenome
		sys.stderr.write(".")
	del genesAnc.dicGenes
	
	# Liste des positions des genes ancestraux dans les genomes modernes
	print >> sys.stderr, " Extraction des positions ...",
	lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
	for e in fils:
		locations[e] = [[] for x in lstGenesAnc]
	for ianc in xrange(len(lstGenesAnc)):
		for g in lstGenesAnc[ianc].names:
			if g not in geneBank.dicGenes:
				continue
			(e,c,i) = geneBank.dicGenes[g]
			locations[e][ianc].append( (c,i,geneBank.dicEspeces[e].lstGenes[c][i].strand) )
	print >> sys.stderr, "OK"


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
			
			utils.myDiags.addDiag(diagEntry[tmp], dd, [(e1,c1), (e2,c2)] )
			
	
	#def addDiag(repos, diag, appar ):
	#
	#	diags = repos[0]
	#	dic = repos[1]
	#	lst = set(utils.myMaths.flatten([dic[x] for x in diag if x in dic]))
	#	flag = False
	#	dd = diag[:]
	#	dd.reverse()
	#	for j in lst:
	#		if utils.myMaths.issublist(diag, diags[j][0]) or utils.myMaths.issublist(dd, diags[j][0]):
	#			diags[j][1].update(appar)
	#			flag = True
	#		elif utils.myMaths.issublist(diags[j][0], diag) or utils.myMaths.issublist(diags[j][0], dd):
	#			diags[j] = ([], set([]))
	#	if not flag:
	#		n = len(diags)
	#		diags.append( (diag,set(appar)) )
	#		for x in diag:
	#			if x not in dic:
	#				dic[x] = []
	#			dic[x].append(n)

	
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

def cutNodes(diags):
	voisins = {}
	for (d,_) in diags:
		if len(d) == 1 and d[0] not in voisins:
			voisins[d[0]] = []
			continue
		for i in xrange(1,len(d)):
			x = d[i-1]
			y = d[i]
			voisins[x] = voisins.get(x,[]) + [y]
			voisins[y] = voisins.get(y,[]) + [x]
	
	for x in voisins:
		voisins[x] = len(set(voisins[x]))

	res = []
	for (d,orig) in diags:
		curr = []
		while len(d) > 0:
			x = d.pop(0)
			if voisins[x] < 3:
				curr.append(x)
			else:
				if len(curr) >= 2:
					res.append( (curr,orig.copy()) )
				curr = []
		if len(curr) >= 2:
			res.append( (curr,orig.copy()) )
	return res


def combinDiags(anc, diags):
	combin = utils.myTools.myCombinator([])
	fils = phylTree.getSpecies(anc)
	for (i,j) in utils.myTools.myMatrixIterator(len(diags), len(diags), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		combin.addLink([i])
		(_,orig1) = diags[i]
		(_,orig2) = diags[j]
		commun = [x for x in orig1.intersection(orig2)]
		filsOK = 0
		outgroupOK = 0
		for (k,l) in utils.myTools.myMatrixIterator(len(commun), len(commun), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			(e1,c1) = commun[k]
			(e2,c2) = commun[l]
			if phylTree.getFirstParent(e1,e2) == anc:
				filsOK += 1
			if ((e1 in fils) and (e2 not in fils)) or ((e2 in fils) and (e1 not in fils)):
				outgroupOK += 1
		if (outgroupOK > 0) and (filsOK > 0):
			combin.addLink([i,j])
	res = []
	for g in combin:
		diag = utils.myMaths.unique(utils.myMaths.flatten([diags[i][0] for i in g]))
		orig = utils.myMaths.unique(utils.myMaths.flatten([diags[i][1] for i in g]))
		res.append( (diag,orig) )
	return res


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("cutNodes",bool,False), ("combinSameChr",bool,False), \
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
diagEntry = dict( [(anc, ([],{})) for anc in phylTree.items] )
calcDiags()

# Plus besoin de ca ...
del genomes
del locations

for anc in diagEntry:
	
	print >> sys.stderr, "Traitement de %s ..." % anc,
	lst = diagEntry[anc][0]
	
	if options["cutNodes"]:
		print >> sys.stderr, "Coupure sur les noeuds ...",
		lst = cutNodes(lst)
	
	if options["combinSameChr"]:
		print >> sys.stderr, "Combinaisons ...",
		lst = combinDiags(anc, lst)
	
	for (d,l) in lst:
		if len(d) == 0:
			continue
		print anc, " ".join([str(x) for x in d])
		#print "%s\t%s\t%s" % (anc, " ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l]))
	
	print >> sys.stderr, "OK"

