#! /users/ldog/muffato/python -OO

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
	if useOutgroups == 2:
		x = [d for (e,d) in distEsp.iteritems() if e in ancSpecies]
		if len(x) > 0:
			x = max(x)
			for e in ancOutgroupSpecies:
				if (e in distEsp) and (distEsp[e] > x):
					del distEsp[e]

	# On calcule par une moyenne les autres distances
	return calcDist(distEsp)

def loadDiagsFile(name):
	f = utils.myTools.myOpenFile(name, 'r')
	diags = utils.myTools.defaultdict(list)
	for l in f:
		t = l.split('\t')
		diags[utils.myGenomes.commonChrName(t[0])].append( ([int(x) for x in t[1].split()],[int(x) for x in t[2].split()]) )
	return diags

#
# Reecrit le genome a trier comme suite des positions des genes dans les genomes modernes
#
def rewriteGenome():

	# On etend la liste des genes ancestraux pour utiliser les outgroup en remontant l'arbre jusqu'a la racine
	dicOutgroupGenes = utils.myTools.defaultdict(set)
	if useOutgroups > 0:
		anc = options["ancestr"]
		while anc in phylTree.parent:
			speciesAlreadyUsed = frozenset(phylTree.species[anc])
			(anc,_) = phylTree.parent[anc]
			# Le genome de l'ancetre superieur
			tmpGenesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
			del tmpGenesAnc.dicGenes
			# Chaque gene
			for g in tmpGenesAnc:
				# Les positions dans les genomes qu'on a charge (on evite le nom FAMxx)
				newGenes = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
				# On se restreint aux outgroup
				newGenes = [x for x in newGenes if x[0] not in speciesAlreadyUsed]
				# On enregistre le lien entre les genes du genome a ordonner et les genes des outgroups
				for x in genesAnc.getPosition(g.names):
					dicOutgroupGenes[x].update(newGenes)
			del tmpGenesAnc


	# On reecrit le genome a trier
	# Chaque gene devient la liste de ses positions dans chaque espece
	genome = {}
	for (c,dd) in diags.iteritems():
		genome[c] = []
		for (d,_) in dd:
			d2 = []
			for i in d:
				tmp = [phylTree.dicGenes[s] for s in genesAnc.lstGenes[None][i].names if s in phylTree.dicGenes]
				tmp.extend(dicOutgroupGenes.get( (None,i), []))
				d2.append(tmp)
			genome[c].append(d2)
	del phylTree.dicGenes
	return genome



# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestralDiags", "phylTree.conf"], \
	[("seuilMaxDistInterGenes",int,100000), ("nbDecimales",int,2), ("infiniteDist",int,1000000), ("notConstraintPenalty",float,10000), \
	("ancestr",str,""), ("nbConcorde",int,1), ("useOutgroups",int,[0,1,2]), ("withConcordeOutput",bool,False),\
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
diags = loadDiagsFile(noms_fichiers["genomeAncestralDiags"])
genesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[options["ancestr"]])
nbConcorde = max(1, options["nbConcorde"])
mult = pow(10, options["nbDecimales"])
seuil = options["seuilMaxDistInterGenes"]
pen = str(options["infiniteDist"] * mult)
add = options["notConstraintPenalty"]
anc = options["ancestr"]
ancSpecies = frozenset(phylTree.species[anc])
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]

phylTree.initCalcDist(anc, useOutgroups != 0)
calcDist = phylTree.calcDist
newGenome = rewriteGenome()
concordeInstance = utils.concorde.ConcordeLauncher()

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for (c,tab) in newGenome.iteritems():

	# Ecrit la matrice du chromosome c et lance concorde
	n = len(tab)
	print >> sys.stderr, "\n- Chromosome %s (%d diagonales) -" % (c, n)
	
	def f(i1, i2):

		#print
		#print
		#print
		#print
		#print
		#print "*** NEW CAL (%d/%d) ***" % (i1,i2)
		#print

		d1 = i1/2
		d2 = i2/2
		if d1 == d2:
			return 0
		
		diag1 = tab[d1]
		if i1 != 2*d1:
			diag1.reverse()
		diag2 = tab[d2]
		if i2 != 2*d2:
			diag2.reverse()

		#print "LENGTHS", len(diag1), len(diag2)
		#print "DIAG1", diag1
		#print "DIAG2", diag2
		#print
		dist = []

		for (i1,g1) in enumerate(diag1):
			for (i2,g2) in enumerate(diag2):
				#print "POSITIONS", i1, i2
				#print "G1", g1
				#print "G2", g2
				y = distInterGenes(g1, g2)
				#print "DIST", y
				#print
				if y != None:
					#print "CORRECTION",  y - i1 - i2,
					y = y - i1 - i2
					if y == 0:
						#print "FINAL"
						return int(mult*add)
					elif y >= 0:
						#print "OK"
						dist.append(y)
					#else:
					#	print "NO"

		if len(dist) == 0:
			#print "FINAL NO"
			return pen
		else:
			#print "FINAL", min(dist)
			return int(mult*(min(dist)+add))
	
	lstTot = concordeInstance.doConcorde(2*n, f, nbConcorde, options["withConcordeOutput"])
	lstTot = utils.myMaths.unique(lstTot)
	if len(lstTot) == 0:
		lstTot = [range(2*n)]
	res = lstTot[0][:]
	tab = diags[c]
		
	while len(res) > 0:
		i1 = res.pop(0)
		i2 = res.pop(0)
		
		diag = tab[i1/2][0]
		strand = tab[i1/2][1]

		if i1/2 != i2/2:
			print >> sys.stderr, "!"
		
		print utils.myTools.printLine([c, 1-2*int(i1>i2), utils.myTools.printLine(diag, " "), utils.myTools.printLine(strand, " ")])
		continue

		#el
		if i1 > i2:
			# La diagonale est a inverser
			diag.reverse()
			strand.reverse()
			for i in xrange(len(strand)):
				strand[i] *= -1
		# On affiche les genes
		for (i,g) in enumerate(diag):
			print c, strand[i], " ".join(genesAnc.lstGenes[None][g].names)
		
	#for sol in lstTot:
	#	print "# .%s" % c, utils.myTools.printLine(sol, " ")
	
	print >> sys.stderr, len(lstTot), "solutions OK"

