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
		global diagEntry

		if len(d1) < options["minimalLength"]:
			return
		
		dn1 = [g1.lstGenes[c1][i].names[0] for i in d1]
		dn2 = [g2.lstGenes[c2][i].names[0] for i in d2]

		for tmp in toStudy:
			diagEntry[tmp].append( ((e1,c1,dn1), (e2,c2,dn2)) )
	
	# Ecrire un genome en suite de genes ancestraux
	def translateGenome(genome):
		newGenome = {}
		for c in genome.lstChr:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
			if options["keepOnlyOrthos"]:
				newGenome[c] = [x for x in newGenome[c] if x[0] != -1]
		return newGenome

	# Les noeuds de l'arbre entre l'ancetre et les especes actuelles
	toStudy = [phylTree.getFirstParent(e1,e2)]
	for tmp in phylTree.items:
		s = phylTree.species[tmp]
		if (e1 in s) and (e2 not in s):
			toStudy.append( tmp )
		elif (e2 in s) and (e1 not in s):
			toStudy.append( tmp )

	# Chargement des orthologues
	#calcDiags(phylTree.commonNames[phylTree.listSpecies[i]][0], phylTree.commonNames[phylTree.listSpecies[j]][0])
	genesAnc = utils.myGenomes.GenomeFromOrthosList(options["orthosFile"] % (phylTree.commonNames[e1][0],phylTree.commonNames[e2][0]))
	newLoc = [[] for x in xrange(len(genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]))]
	del genesAnc.lstGenes

	g1 = phylTree.dicGenomes[e1]
	g2 = phylTree.dicGenomes[e2]
	newGen = translateGenome(g1)
	tmp = translateGenome(g2)
	
	for c in g2.lstChr:
		for i in xrange(len(tmp[c])):
			(ianc,s) = tmp[c][i]
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )
	
	utils.myDiags.iterateDiags(newGen, newLoc, options["fusionThreshold"], options["sameStrand"], combinDiag)


# Cherche si la diagonale est localisee sur un chromosome d'une espece, non ordonnee
def findNewSpecies(d, esp, anc):
	
	# Les especes pour lesquelles la diagonale n'a pas ete detectee
	lstEsp = set(phylTree.listSpecies).difference([e for (e,_) in esp])
	res = []
	for e in lstEsp:
		# L'ancetre commun
		a = phylTree.getFirstParent(anc, e)
		# On va intersecter les chromosomes des orthologues de chaque gene
		poss = set([])
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
	("extractLongestPath",bool,True), ("searchUndetectedSpecies",bool,True), \
	("genesFile",str,"~/work/data/genes/full/genes.%s.list.bz2"), \
	("orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince("Euteleostomi", options["genesFile"])

# La structure qui accueillera les diagonales
diagEntry = dict( [(anc, []) for anc in phylTree.items] )
# On compare toutes les especes entre elles
for (i,j) in utils.myTools.myMatrixIterator(len(phylTree.listSpecies), len(phylTree.listSpecies), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	print >> sys.stderr, '+',
	calcDiags(phylTree.listSpecies[i], phylTree.listSpecies[j])

# On a besoin des genes ancestraux
genesAnc = {}
for anc in diagEntry:
	genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % (anc.replace('/', '_').replace(' ', '_')))

# Traitement final
for anc in diagEntry:
	
	s = []
	if options["extractLongestPath"]:
	
		print >> sys.stderr, "Extraction des chevauchements les plus longs de %s " % anc,
		lst = utils.myDiags.extractLongestOverlappingDiags(diagEntry[anc], genesAnc[anc])

		print >> sys.stderr, "OK (%d -> %d) ... Impression ..." % (len(diagEntry[anc]), len(lst)),
		for (l,d,esp) in lst:
			s.append( l )
			supp = findNewSpecies(d, esp, anc)
			print "%s\t%d\t%s\t%s\t%s" % \
			(anc, l, " ".join([str(x) for x in d]), "|".join(["%s/%s" % (e,c) for (e,c) in esp]), "|".join(["%s/%s" % (e,c) for (e,c) in supp]))
	
	else:
	
		lst = diagEntry[anc]
		print >> sys.stderr, "Impression des %d diagonales de %s ..." % (len(lst),anc),
		for ((e1,c1,d1),(e2,c2,d2)) in lst:
			s.append( len(d1) )
			print "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (anc, len(d1), e1,c1," ".join(d1), e2,c2," ".join(d2))
		
	ss = sum(s)
	if len(lst) == 0:
		print >> sys.stderr, ss, "%.2f" % 0, 0, "OK"
		continue
	print >> sys.stderr, ss, "%.2f" % (float(ss)/float(len(lst))), max(s), "OK"



