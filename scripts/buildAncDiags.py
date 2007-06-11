#! /users/ldog/muffato/python -OO

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
import utils.myBioObjects
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags


#############
# FONCTIONS #
#############


def getLongestDiags(oldDiags):

	dic = {}
	diags = []
	combin = utils.myTools.myCombinator([])
	for i in xrange(len(oldDiags)):
		((_,_,d1),(_,_,d2),_) = oldDiags[i]
		((e1,c1,d1),(e2,c2,d2),da) = oldDiags[i]
		da1 = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d1]
		if "" in da1:
			diags.append( [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d2] )
		else:
			diags.append(da1)
		if "" in diags[-1]:
			print >> sys.stderr, e1,c1
			print >> sys.stderr, e2,c2
			print >> sys.stderr, d1
			print >> sys.stderr, d2
			print >> sys.stderr, da
		for s in d1+d2:
			if s not in dic:
				dic[s] = []
			dic[s].append(i)
		combin.addLink([i])
	
	for s in dic:
		combin.addLink(dic[s])
	
	newDiags = []
	for g in combin:
		if options["extractLongestPath"]:
			gr = utils.myDiags.DiagGraph([diags[i] for i in g])
			gr.reduceGraph()
			#print >> sys.stderr, "%d/%d->%d" % (len(gr.sommets),len(g),len(gr.newSommets)),
			#gr.printGraph()
			#gr.printReducedGraph()
		else:
			gr = utils.myDiags.WeightedDiagGraph([diags[i] for i in g])
		for res in gr.getBestDiags():
			if len(res) < options["minimalLength"]:
				continue
			ok = set()
			for i in g:
				d = diags[i]
				for j in xrange(len(d)-1):
					try:
						i1 = res.index(d[j])
						i2 = res.index(d[j+1])
					except ValueError:
						continue
					if abs(i1-i2) != 1:
						continue
					ok.add( (oldDiags[i][0][0],oldDiags[i][0][1]) )
					ok.add( (oldDiags[i][1][0],oldDiags[i][1][1]) )
					break
			newDiags.append( (len(res), res, list(ok)) )
	return newDiags
	



# Cherche si la diagonale est localisee sur un chromosome d'une espece, non ordonnee
def findNewSpecies(d, esp, anc):
	
	# Les especes pour lesquelles la diagonale n'a pas ete detectee
	lstEsp = set(listSpecies).difference([e for (e,_) in esp])
	res = []
	for e in lstEsp:
		# L'ancetre commun
		a = phylTree.dicParents[anc][e]
		# On va intersecter les chromosomes des orthologues de chaque gene
		poss = set()
		for i in d:
			
			# Le gene chez l'ancetre commun
			g = genesAnc[a].getPosition(genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr][i])
			# Cas d'un gene specifique de la lignee
			if len(g) == 0:
				continue
			# Le gene dans l'autre espece
			tmp = [c for (c,_) in dicGenomes[e].getPosition(genesAnc[a].lstGenes[g[0][0]][g[0][1]])]
			# Gene non trouve, on passe au suivant
			if len(tmp) == 0:
				continue
			# Sinon, on intersecte
			if len(poss) == 0:
				poss.update(tmp)
			else:
				poss.intersection_update(tmp)
				# Utile de continuer ?
				if len(poss) == 0:
					break
		# S'il en reste un, c'est gagne !
		if len(poss) != 0:
			res.append( (e,poss.pop()) )
	
	return res


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOnlyOrthos",bool,False),
	("useOutgroups",bool,False), ("target",str,""), \
	("showProjected",bool,False), ("showAncestral",bool,True), ("searchUndetectedSpecies",bool,True), \
	("extractLongestPath",bool,False), ("cutLongestPath",bool,False), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("genesAncFile",str,"~/work/ancestralGenomes/Genome.%s.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Les especes a utiliser
dicGenomes = {}
tmp = options["target"].split(',')
if len(options["target"]) == 0:
	print >> sys.stderr, "Aucune cible indiquee pour l'extraction des diagonales"
	sys.exit(1)
else:
	listSpecies = []
	for x in tmp:
		if x[0] != '.':
			listSpecies.extend(phylTree.species[x])
			for e in phylTree.species[x]:
				dicGenomes[e] = utils.myGenomes.EnsemblGenome(options["genesFile"] % phylTree.fileName[e])
		else:
			listSpecies.append(x[1:])
			dicGenomes[x[1:]] = utils.myGenomes.AncestralGenome(options["genesAncFile"] % phylTree.fileName[x[1:]], chromPresents=True)

# Les outgroup du noeud le plus ancien
if options.useOutgroups:
	target = listSpecies[0]
	for e in listSpecies:
		target = phylTree.dicParents[target][e]
	listSpecies += phylTree.outgroupSpecies[target]
	for e in phylTree.outgroupSpecies[target]:
		dicGenomes[e] = utils.myGenomes.EnsemblGenome(options["genesFile"] % phylTree.fileName[e])

# La liste des ancetres edites
tmp = set()
for (e1,e2) in utils.myTools.myMatrixIterator(listSpecies, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
	tmp.update(phylTree.dicLinks[e1][e2][1:-1])
	tmp.add(phylTree.dicParents[e1][e2])
diagEntry = {}
genesAnc = {}
for anc in tmp:
	diagEntry[anc] = []
	genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[anc])

	
# La fonction qui permet de stocker les diagonales sur les ancetres
def storeDiag(da):
	for anc in toStudy:
		diagEntry[anc].append( da )
	
# On compare toutes les especes entre elles
for (e1,e2) in utils.myTools.myMatrixIterator(listSpecies, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
	toStudy = set(phylTree.dicLinks[e1][e2][1:-1] + [phylTree.dicParents[e1][e2]])
	utils.myDiags.calcDiags(e1, e2, dicGenomes[e1], dicGenomes[e2], genesAnc[phylTree.dicParents[e1][e2]], storeDiag, options["minimalLength"], \
		options["fusionThreshold"], options["sameStrand"] and (e1 not in genesAnc) and (e2 not in genesAnc), options["keepOnlyOrthos"])

# Traitement final
for anc in diagEntry:

	if options["showProjected"]:
		lst = diagEntry[anc]
		print >> sys.stderr, "Impression des %d diagonales projetees de %s ..." % (len(lst),anc),
		s = []
		for ((e1,c1,d1),(e2,c2,d2),_) in lst:
			s.append( len(d1) )
			print '\t'.join([anc, str(len(d1)), e1,str(c1)," ".join(d1), e2,str(c2)," ".join(d2)])
		print >> sys.stderr, utils.myMaths.myStats(s), "OK"


	if options["showAncestral"]:
	
		if options["extractLongestPath"] or options["cutLongestPath"]:
			print >> sys.stderr, "Extraction des chevauchements les plus longs de %s ..." % anc,
			lst = getLongestDiags(diagEntry[anc])
			print >> sys.stderr, "OK (%d -> %d)" % (len(diagEntry[anc]), len(lst))
		else:
			lst = []
			for ((e1,c1,d1),(e2,c2,d2),_) in diagEntry[anc]:
				da = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d1]
				if "" in da:
					da = [genesAnc[anc].dicGenes[s][1] for s in d2]
				lst.append( (len(da), da, [(e1,c1),(e2,c2)]) )

		
		print >> sys.stderr, "Impression des %d diagonales ancestrales de %s ..." % (len(lst),anc),
		s = []
		for (l,d,esp) in lst:
			s.append( l )
			ss = '\t'.join([anc, str(l), " ".join([str(x) for x in d]), "|".join(["%s/%s" % (e,c) for (e,c) in esp])])
			if options["searchUndetectedSpecies"]:
				supp = findNewSpecies(d, esp, anc)
				ss += '\t' + '|'.join(["%s/%s" % (e,c) for (e,c) in supp])
			print ss
	
		print >> sys.stderr, utils.myMaths.myStats(s), "OK"

