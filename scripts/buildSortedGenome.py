#! /users/ldog/muffato/python

__doc__ = """
Trie les gens selon l'ordre consensus issu des especes de l'arbre phylogenetique
"""

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
	if useOutgroups == "onlyIfBetter":
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
	if useOutgroups != "no":
		anc = arguments["ancestr"]
		while anc in phylTree.parent:
			(anc,_) = phylTree.parent[anc]
			# Le genome de l'ancetre superieur
			tmpGenesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
			del tmpGenesAnc.dicGenes
			# Chaque gene
			for g in tmpGenesAnc:
				# Les positions dans les genomes qu'on a charge (on evite le nom FAMxx)
				newGenes = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
				# On se restreint aux outgroup
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
arguments = utils.myTools.checkArgs( \
	[("genomeAncestral",file), ("phylTree.conf",file)], \
	[("ancestr",str,""), \
	("seuilMaxDistInterGenes",int,1000000), ("nbDecimales",int,2), ("infiniteDist",int,1000000), ("notConstraintPenalty",float,0), \
	("useOutgroups",str,["no","always","onlyIfBetter"]), ("newParsimonyScoring",bool,False), \
	("nbConcorde",int,1), ("withConcordeOutput",bool,False), ("searchAllSol",bool,False), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogentique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
if arguments["ancestr"] not in phylTree.allNames:
	print >> sys.stderr, "Can't retrieve the order of -%s- " % arguments["ancestr"]
	sys.exit(1)
# On charge les genomes
useOutgroups = arguments["useOutgroups"]
if useOutgroups != "no":
	root = None
else:
	root = arguments["ancestr"]
phylTree.loadAllSpeciesSince(root, arguments["genesFile"], storeGenomes = False)
genesAnc = utils.myGenomes.Genome(arguments["genomeAncestral"])
nbConcorde = max(1, arguments["nbConcorde"])
mult = pow(10, arguments["nbDecimales"])
seuil = arguments["seuilMaxDistInterGenes"]
pen = str(int(mult * arguments["infiniteDist"]))
add = arguments["notConstraintPenalty"]
anc = arguments["ancestr"]
ancSpecies = phylTree.species[anc]
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]

phylTree.initCalcDist(anc, useOutgroups != "no")
if arguments["newParsimonyScoring"]:
	calcDist = lambda distEsp: phylTree.calcWeightedValue(distEsp, -1, anc)[2]
else:
	calcDist = phylTree.calcDist
newGenome = rewriteGenome()

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for (c,tab) in newGenome.iteritems():

	# Ecrit la matrice du chromosome c et lance concorde
	n = len(tab)
	print >> sys.stderr, "\n- Chromosome %s (%d genes) -" % (c, n)
	
	# La matrice des distances intergenes
	mat = [[None] * n for i in xrange(n)]
	
	def f(i1, i2):
		mat[i1][i2] = mat[i2][i1] = y = distInterGenes(tab[i1], tab[i2])
		# L'astuce est que None < 0
		if y < 0:
			return pen
		elif y == 1:
			return 0
		else:
			return int(mult*y+add)

	lstTot = utils.concorde.doConcorde(n, f, nbConcorde, arguments["withConcordeOutput"])
	lstTot = utils.myMaths.unique(lstTot)
	if len(lstTot) == 0:
		lstTot = [range(n)]

	for (i,g) in enumerate(lstTot[0]):
		print c,
		if nbConcorde > 1:
			print len(set([s[i] for s in lstTot])),
		print " ".join(genesAnc.lstGenes[c][g].names)
	
	if nbConcorde > 1:
		for sol in solUniq:
			print "# .", c, utils.myTools.printLine(sol)
	print >> sys.stderr, len(solUniq), "solutions"

	if arguments["searchAllSol"]:
		print >> sys.stderr, "Etude des solutions alternatives ...",
		nb = 0
		res = lstTot[0]
		tmp = [(i,res[i],res[i+1],mat[res[i]][res[i+1]]) for i in xrange(n-1)]
		for (i,xi,xip,diip) in tmp:
			for (j,xj,xjp,djjp) in tmp[i+2:]:
				dij = mat[xi][xj]
				if dij < 0:
					continue
				dijpp = mat[xip][xjp]
				if dijpp < 0:
					continue
				if int(mult*(diip+djjp)) == int(mult*(dij+dijpp)):
					# On a detecte une inversion potentielle
					print "# ... %d / %d ... %d / %d ..." % (i,i+1,j,j+1)
					nb += 1
		
		print >> sys.stderr, nb, "inversions possibles"
	
