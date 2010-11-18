#!/usr/bin/env python2

__doc__ = """
Trie les gens selon l'ordre consensus issu des especes de l'arbre phylogenetique
"""

import sys
import collections

import utils.myFile
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

def loadDiagsFile(name):
	f = utils.myFile.openFile(name, 'r')
	diags = collections.defaultdict(list)
	for l in f:
		t = l.split('\t')
		diags[t[0]].append( ([int(x) for x in t[1].split()],[int(x) for x in t[2].split()]) )
		#diags[utils.myGenomes.commonChrName(t[0])].append( ([int(x) for x in t[1].split()],[int(x) for x in t[2].split()]) )
	return diags

#
# Reecrit le genome a trier comme suite des positions des genes dans les genomes modernes
#
def rewriteGenome():

	# On etend la liste des genes ancestraux pour utiliser les outgroup en remontant l'arbre jusqu'a la racine
	dicOutgroupGenes = collections.defaultdict(set)
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
arguments = utils.myTools.checkArgs( \
	[("genomeAncestralDiags",file), ("phylTree.conf",file), ("ancestr",str)], \
	[("seuilMaxDistInterGenes",int,1000000), ("nbDecimales",int,2), ("infiniteDist",int,1000000), ("notConstraintPenalty",float,0), \
	("useOutgroups",str,["no","always","onlyIfBetter"]), ("newParsimonyScoring",bool,False), \
	("nbConcorde",int,1), ("withConcordeOutput",bool,False), \
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
diags = loadDiagsFile(arguments["genomeAncestralDiags"])
genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[arguments["ancestr"]])
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
	calcDist = lambda distEsp: phylTree.calcDist(distEsp, 0)
newGenome = rewriteGenome()

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for (c,tab) in newGenome.iteritems():

	# Ecrit la matrice du chromosome c et lance concorde
	n = len(tab)
	print >> sys.stderr, "\n- Chromosome %s (%d diagonales) -" % (c, n)
	
	def f(i1, i2):

		d1 = i1/2
		d2 = i2/2
		if d1 == d2:
			return 0
		
		diag1 = tab[d1]
		if i1 != 2*d1:
			diag1 = list(reversed(diag1))
		diag2 = tab[d2]
		if i2 != 2*d2:
			diag2 = list(reversed(diag2))

		dist = []
		for (i1,g1) in enumerate(diag1):
			for (i2,g2) in enumerate(diag2):
				y = distInterGenes(g1, g2)
				if y != None:
					y = y - i1 - i2
					if y == 0:
						return int(mult*add)
					elif y > 0:
						dist.append(y)

		if len(dist) == 0:
			print len(diag1)
			print diag1[0]
			print diag1[-1]
			print len(diag2)
			print diag2[0]
			print diag2[-1]
			print pen
			print
			return pen
		else:
			y = min(dist)
			print len(diag1)
			print diag1[0]
			print diag1[-1]
			print len(diag2)
			print diag2[0]
			print diag2[-1]
			print int(mult*y+add)
			print
			return int(mult*y+add)

	lstTot = utils.concorde.doConcorde(2*n, f, nbConcorde, arguments["withConcordeOutput"])
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
		
		print utils.myFile.myTSV.printLine([c, 1-2*int(i1>i2), utils.myTools.printLine(diag, " "), utils.myTools.printLine(strand, " ")])
		continue

		if i1 > i2:
			# La diagonale est a inverser
			diag.reverse()
			strand.reverse()
			for i in xrange(len(strand)):
				strand[i] *= -1
		# On affiche les genes
		for (i,g) in enumerate(diag):
			print c, strand[i], " ".join(genesAnc.lstGenes[None][g].names)
	
	print >> sys.stderr, len(lstTot), "solutions OK"
