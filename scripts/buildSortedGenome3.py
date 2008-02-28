#! /users/ldog/muffato/python -OO

__doc__ = """
Trie les gens selon l'ordre consensus issu des especes de l'arbre phylogenetique
"""

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.concorde

itStrictUpperMatrix = utils.myTools.myIterator.tupleOnStrictUpperList

#
# Calcule la distance inter-genes moyenne entre les deux genes
#   Utilise la fonction d'inference de valeur de phylTree
#
def distInterGenes(tg1, tg2):

	# Les distances chez chaque espece
	distEsp = {}
	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				# Au dela d'un certain seuil, on considere que l'information n'est plus valable
				# Si les deux genes ont la meme coordonnee, ce sont des duplicats et on n'a pas de distance a comptabiliser
				if (x > 0) and (x <= seuil):
					# On garde la plus petite distance trouvee
					if (e1 not in distEsp) or (x < distEsp[e1]):
						distEsp[e1] = x
	
	# On fait la liste des especes qui presentent une distance de 1
	# On met les 1 dans les noeuds de l'arbre entre ces especes
	lst1Anc = set()
	for (e1,e2) in itStrictUpperMatrix([e for (e,d) in distEsp.iteritems() if d == 1]):
		lst1Anc.update(phylTree.dicLinks[e1][e2])
		if anc in lst1Anc:
			return 1
	for e in lst1Anc:
		distEsp[e] = 1

	# En mode outgroups/2 les outgroups qui montrent une distance plus grande que les fils sont supprimes
	if useOutgroups == 2:
		x = [d for (e,d) in distEsp.iteritems() if e in ancSpecies]
		if len(x) > 0:
			x = max(x)
			for e in ancOutgroupSpecies:
				if (e in distEsp) and (distEsp[e] > x):
					del distEsp[e]

	# On calcule par une moyenne les autres distances
	return calcDist(distEsp)


#
# Reecrit le genome a trier comme suite des positions des genes dans les genomes modernes
#
def rewriteGenome():

	# On etend la liste des genes ancestraux pour utiliser les outgroup en remontant l'arbre jusqu'a la racine
	dicOutgroupGenes = utils.myTools.defaultdict(set)
	if useOutgroups > 0:
		anc = options["ancestr"]
		while anc in phylTree.parent:
			(anc,_) = phylTree.parent[anc]
			# Le genome de l'ancetre superieur
			tmpGenesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
			del tmpGenesAnc.dicGenes
			# Chaque gene
			for g in tmpGenesAnc:
				# Les autres sont les orthologues dans les autres especes
				newGenes = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
				newGenes = [x for x in newGenes if x[0] not in phylTree.species[anc]]
				# On enregistre le lien entre les genes du genome a ordonner et les genes des outgroups
				for x in genesAnc.getPosition(g.names):
					dicOutgroupGenes[x].update(newGenes)
			del tmpGenesAnc


	# On reecrit le genome a trier
	# Chaque gene devient la liste de ses positions dans chaque espece
	genome = {}
	for c in genesAnc.lstChr:
		genome[c] = []
		for (i,g) in enumerate(genesAnc.lstGenes[c]):
			tmp = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
			tmp.extend(dicOutgroupGenes.get( (c,i), []))
			genome[c].append(tmp)
	del phylTree.dicGenes
	return genome



# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestral", "phylTree.conf"], \
	[("ancestr",str,""), ("seuilMaxDistInterGenes",int,1000000), ("nbDecimales",int,2), ("infiniteDist",int,1000000), ("notConstraintPenalty",float,0), \
	("useOutgroups",int,[0,1,2]), ("newParsimonyScoring",bool,False), ("nbConcorde",int,1), ("withConcordeOutput",bool,False), ("searchAllSol",bool,False), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogentique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["ancestr"] not in phylTree.allNames:
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["ancestr"]
	sys.exit(1)
# On charge les genomes
useOutgroups = options["useOutgroups"]
if useOutgroups > 0:
	root = None
else:
	root = options["ancestr"]
phylTree.loadAllSpeciesSince(root, options["genesFile"], storeGenomes = False)
genesAnc = utils.myGenomes.Genome(noms_fichiers["genomeAncestral"])
nbConcorde = max(1, options["nbConcorde"])
mult = pow(10, options["nbDecimales"])
seuil = options["seuilMaxDistInterGenes"]
pen = str(int(mult*options["infiniteDist"]))
add = options["notConstraintPenalty"]
anc = options["ancestr"]
ancSpecies = phylTree.species[anc]
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]

phylTree.initCalcDist(anc, useOutgroups != 0)
if options["newParsimonyScoring"]:
	calcDist = lambda distEsp: phylTree.calcWeightedValue(distEsp, -1, anc)[2]
else:
	calcDist = phylTree.calcDist
genome = rewriteGenome()
concordeInstance = utils.concorde.ConcordeLauncher()

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:

	# Ecrit la matrice du chromosome c et lance concorde
	tab = genome[c]
	n = len(tab)
	print >> sys.stderr, "\n- Chromosome %s (%d genes) -" % (c, n)
	
	# La matrice des distances intergenes
	mat = [[None] * (n+1) for i in xrange(n+1)]
	
	def f(i, j):
		y = distInterGenes(tab[i], tab[j])
		mat[i+1][j+1] = mat[j+1][i+1] = y
		if y < 0:
			y = None
		if y == None:
			return pen
		elif y == 1:
			return 0
		else:
			return int(mult*y+add)

	lstTot = concordeInstance.doConcorde(n, f, nbConcorde, options["withConcordeOutput"])

	if len(lstTot) == 0:
		lstTot = [range(n)]
	for (i,g) in enumerate(lstTot[0]):
		print c,
		if nbConcorde > 1:
			print len(set([s[i] for s in lstTot])),
		print " ".join(genesAnc.lstGenes[c][g].names)
	
	solUniq = utils.myMaths.unique(lstTot)
	if nbConcorde > 1:
		for sol in solUniq:
			print ".", c, utils.myTools.printLine(sol)
	print >> sys.stderr, len(solUniq), "solutions"

	if options["searchAllSol"]:
		print >> sys.stderr, "Etude des solutions alternatives ...",
		nb = 0
		res = lstTot[0]
		tmp = [(i,res[i],res[i+1],mat[res[i]][res[i+1]]) for i in xrange(n)]
		for (i,xi,xip,diip) in tmp:
			for (j,xj,xjp,djjp) in tmp[:max(i-1,0)]:
				dij = mat[xi][xj]
				if dij == None:
					continue
				dijpp = mat[xip][xjp]
				if dijpp == None:
					continue
				if int(mult * (diip+djjp-dij-dijpp)) == 0:
					# On a detecte une inversion potentielle
					print "# ... %d / %d ... %d / %d ..." % (j,j+1,i,i+1)
					nb += 1
		
		print >> sys.stderr, nb, "inversions possibles"
	
