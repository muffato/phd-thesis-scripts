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

def calcDiags(e1, e2):

	
	# La fonction qui permet de stocker les diagonales sur les ancetres
	def combinDiag(c1, c2, d1, d2):
		global diagEntry, statsDiags

		if len(d1) < options["minimalLength"]:
			return
		
		statsDiags.append(len(d1))
		dn1 = [g1.lstGenes[c1][trans1[(c1,i)]].names[0] for i in d1]
		dn2 = [g2.lstGenes[c2][trans2[(c2,i)]].names[0] for i in d2]

		for tmp in toStudy:
			diagEntry[tmp].append( ((e1,c1,dn1), (e2,c2,dn2)) )
	
	# Ecrire un genome en suite de genes ancestraux
	def translateGenome(genome):
		newGenome = {}
		transNewOld = {}
		for c in genome.lstChr + genome.lstScaff:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
			if options["keepOnlyOrthos"]:
				tmp = [x for x in newGenome[c] if x[0] != -1]
			else:
				tmp = newGenome[c]
			last = 0
			for i in xrange(len(tmp)):
				x = tmp[i]
				new = newGenome[c].index(x, last)
				transNewOld[(c,i)] = new
				last = new + 1
			newGenome[c] = tmp
					
		return (newGenome,transNewOld)


	# Chargement des orthologues
	genesAnc = utils.myGenomes.EnsemblOrthosListGenome(options["orthosFile"] % (phylTree.fileName[e1],phylTree.fileName[e2]), ancFilter=[phylTree.getFirstParent(e1,e2)])
	newLoc = [[] for x in xrange(len(genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]))]
	del genesAnc.lstGenes
	
	# Les noeuds de l'arbre entre l'ancetre et les especes actuelles
	toStudy = utils.myMaths.flatten(phylTree.dicLinks[(e1,e2)])
	global diagEntry
	for tmp in toStudy:
		diagEntry[tmp] = diagEntry.get(tmp, [])

	g1 = phylTree.dicGenomes[e1]
	g2 = phylTree.dicGenomes[e2]
	(newGen,trans1) = translateGenome(g1)
	(tmp,trans2) = translateGenome(g2)
	
	for c in g2.lstChr + g2.lstScaff:
		for i in xrange(len(tmp[c])):
			(ianc,s) = tmp[c][i]
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )

	global statsDiags
	print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
	statsDiags = []
	utils.myDiags.iterateDiags(newGen, newLoc, options["fusionThreshold"], options["sameStrand"], combinDiag)
	print >> sys.stderr, utils.myMaths.myStats(statsDiags)

def getLongestDiags(oldDiags):

	dic = {}
	diags = []
	combin = utils.myTools.myCombinator([])
	for i in xrange(len(oldDiags)):
		d1 = oldDiags[i][0][2]
		d2 = oldDiags[i][1][2]
		da1 = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d1]
		da2 = [genesAnc[anc].dicGenes.get(s,("",""))[1] for s in d2]
		if "" in da1:
			diags.append(da2)
		else:
			diags.append(da1)
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
		a = phylTree.getFirstParent(anc, e)
		# On va intersecter les chromosomes des orthologues de chaque gene
		poss = set()
		for i in d:
			
			# Le gene chez l'ancetre commun
			g = genesAnc[a].getPosition(genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr][i])
			# Cas d'un gene specifique de la lignee
			if len(g) == 0:
				continue
			# Le gene dans l'autre espece
			tmp = [c for (c,_) in phylTree.dicGenomes[e].getPosition(genesAnc[a].lstGenes[utils.myGenomes.Genome.defaultChr][g[0][1]])]
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
	("genesFile",str,"~/work/data/genes/full/genes.%s.list.bz2"), \
	("orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Les especes a utiliser
tmp = options["target"].split(',')
if len(tmp) == 0:
	print >> sys.stderr, "Aucune cible indiquee pour l'extraction des diagonales"
	sys.exit(1)
elif len(tmp) == 1:
	listSpecies = phylTree.species[tmp[0]]
else:
	listSpecies = [phylTree.officialName[x] for x in tmp]
if options.useOutgroups:
	listSpecies += phylTree.outgroupSpecies[target]
phylTree.loadSpeciesFromList(listSpecies, options.genesFile)

# On compare toutes les especes entre elles
diagEntry = {}
for (i,j) in utils.myTools.myMatrixIterator(len(listSpecies), len(listSpecies), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	calcDiags(listSpecies[i], listSpecies[j])

# On a besoin des genes ancestraux
if options["showAncestral"]:
	genesAnc = {}
	for anc in diagEntry:
		genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[anc])

# Traitement final
for anc in diagEntry:

	if options["showProjected"]:
		lst = diagEntry[anc]
		print >> sys.stderr, "Impression des %d diagonales projetees de %s ..." % (len(lst),anc),
		s = []
		for ((e1,c1,d1),(e2,c2,d2)) in lst:
			s.append( len(d1) )
			print '\t'.join([anc, str(len(d1)), e1,str(c1)," ".join(d1), e2,str(c2)," ".join(d2)])
		ss = sum(s)
		if len(lst) == 0:
			print >> sys.stderr, ss, "%.2f" % 0, 0, "OK"
			continue
		print >> sys.stderr, ss, "%.2f" % (float(ss)/float(len(lst))), max(s), "OK"


	if options["showAncestral"]:
	
		if options["extractLongestPath"] or options["cutLongestPath"]:
			print >> sys.stderr, "Extraction des chevauchements les plus longs de %s ..." % anc,
			lst = getLongestDiags(diagEntry[anc])
			print >> sys.stderr, "OK (%d -> %d)" % (len(diagEntry[anc]), len(lst))
		else:
			lst = []
			for ((e1,c1,d1),(e2,c2,d2)) in diagEntry[anc]:
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
	
		ss = sum(s)
		if len(lst) == 0:
			print >> sys.stderr, ss, "%.2f" % 0, 0, "OK"
			continue
		print >> sys.stderr, ss, "%.2f" % (float(ss)/float(len(lst))), max(s), "OK"

