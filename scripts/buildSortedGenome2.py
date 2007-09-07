#! /users/ldog/muffato/python -OO

__doc__ = """
Trie les gens selon l'ordre consensus issu des especes de l'arbre phylogenetique
"""

# Librairies
import os
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.concorde

itStrictUpperMatrix = utils.myTools.myIterator.tupleOnStrictUpperList


#############
# FONCTIONS #
#############

#
# Calcule la distance inter-genes moyenne entre les deux genes
#   Utilise la fonction d'inference de valeur de phylTree
#
def distInterGenes(tg1, tg2, seuil):

	# Les distances chez chaque espece
	distEsp = tabDist[:]

	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				# Au dela d'un certain seuil, on considere que l'information n'est plus valable
				if (seuil <= 0) or (x <= seuil):
					# On garde la plus petite distance trouvee
					if distEsp[e1] == None:
						distEsp[e1] = x
					else:
						distEsp[e1] = min(distEsp[e1], x)
	
	# On fait la liste des especes qui presentent une distance de 1
	lst1Esp = [phylTree.allNames[i] for (i,d) in enumerate(distEsp) if d == 1]
	
	# On met les 1 dans les noeuds de l'arbre entre les especes
	lst1Anc = set()
	for (e1,e2) in itStrictUpperMatrix(lst1Esp):
		lst1Anc.update(phylTree.dicLinks[e1][e2])
		if anc in lst1Anc:
			return 1
	for a in lst1Anc:
		distEsp[phylTree.indNames[a]] = 1

	# En mode outgroups/2 les outgroups qui montrent une distance plus grande que les fils sont supprimes
	if options["useOutgroups"] == 2:
		# HACK: max(None,x) renvoie x
		m = max([distEsp[e] for e in ancSpecies])
		if m != None:
			for e in ancOutgroupSpecies:
				if distEsp[e] > m:
					distEsp[e] = None

	# On calcule par une moyenne les autres distances
	return phylTree.calcDist2(distEsp)


#
# Reecrit le genome a trier comme suite des positions des genes dans les genomes modernes
#
def rewriteGenome():

	# On etend la liste des genes ancestraux pour utiliser les outgroup en remontant l'arbre jusqu'a la racine
	dicOutgroupGenes = {}
	if options["useOutgroups"] > 0:
		anc = options["ancestr"]
		while anc in phylTree.parent:
			anc = phylTree.parent[anc]
			# Le genome de l'ancetre superieur
			tmpGenesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
			del tmpGenesAnc.dicGenes
			# Chaque gene
			for g in tmpGenesAnc:
				ianc = set()
				newGenes = []
				for s in g.names:
					# On se sert des genes qu'on retrouve dans le genome a ordonner comme ancre
					if s in genesAnc.dicGenes:
						ianc.add(genesAnc.dicGenes[s])
					# Les autres sont les orthologues dans les autres especes
					elif s in phylTree.dicGenes:
						tt = phylTree.dicGenes[s]
						if tt[0] not in phylTree.species[options["ancestr"]]:
							newGenes.append(tt)
				# On enregistre le lien entre les genes du genome a ordonner et les genes des outgroups
				for i in ianc:
					if i in dicOutgroupGenes:
						dicOutgroupGenes[i].update(newGenes)
					else:
						dicOutgroupGenes[i] = set(newGenes)
			del tmpGenesAnc


	# On reecrit le genome a trier
	# Chaque gene devient la liste de ses positions dans chaque espece
	genome = {}
	for c in genesAnc.lstChr:
		genome[c] = []
		for (i,g) in enumerate(genesAnc.lstGenes[c]):
			tmp = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
			tmp.extend(dicOutgroupGenes.get( (c,i), []))
			genome[c].append( [(phylTree.indNames[e],ch,i) for (e,ch,i) in tmp] )
	del phylTree.dicGenes
	return genome



# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestral", "phylTree.conf"], \
	[("ancestr",str,""), ("seuilMaxDistInterGenes",float,0), ("nbDecimales",int,2), ("penalite",int,1000000), \
	("useOutgroups",int,[0,1,2]), ("nbConcorde",int,1), ("withConcordeOutput",bool,False), ("withConcordeStats",bool,False),\
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogentique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["ancestr"] not in (phylTree.listAncestr + phylTree.listSpecies):
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["ancestr"]
	sys.exit(1)
# On charge les genomes
if options["useOutgroups"] > 0:
	phylTree.loadAllSpeciesSince(None, options["genesFile"])
else:
	phylTree.loadAllSpeciesSince(options["ancestr"], options["genesFile"])
del phylTree.dicGenomes
genesAnc = utils.myGenomes.Genome(noms_fichiers["genomeAncestral"])
nbConcorde = max(1, options["nbConcorde"])
mult = pow(10, options["nbDecimales"])
seuil = options["seuilMaxDistInterGenes"]
pen = str(int(mult*options["penalite"]))
anc = options["ancestr"]
ancSpecies = [phylTree.indNames[e] for e in phylTree.species[anc]]
ancOutgroupSpecies = [phylTree.indNames[e] for e in phylTree.outgroupSpecies[anc]]
tabDist = [None] * len(phylTree.allNames)
nbNames = len(tabDist)

phylTree.initCalcDist2(anc, options["useOutgroups"] != 0)
genome = rewriteGenome()
concordeInstance = utils.concorde.ConcordeLauncher()

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:

	# Ecrit la matrice du chromosome c et lance concorde
	tab = genome[c]
	n = len(tab)
	
	print >> sys.stderr, "\n- Chromosome %s (%d genes) -" % (c, n)
	
	def f(i, j):
		y = distInterGenes(tab[i], tab[j], seuil)
		if y == None:
			return pen
		elif y == 1:
			return 0
		else:
			return int(mult*y)
	
	#lstTot = concordeInstance.doConcorde(n, f, nbConcorde, options["withConcordeOutput"])
	lstTot = concordeInstance.doConcorde(n, f, 0, options["withConcordeOutput"])

	if len(lstTot) == 0:
		for g in genesAnc.lstGenes[c]:
			print c,
			if nbConcorde > 1:
				print 1,
			print " ".join(g.names)
	else:
		for (i,g) in enumerate(lstTot[0].res):
			print c,
			if nbConcorde > 1:
				print len(set([s.res[i] for s in lstTot])),
			print " ".join(genesAnc.lstGenes[c][g-1].names)
	
	solUniq = utils.myMaths.unique([l.res for l in lstTot])

	if options["withConcordeStats"]:
		for sol in solUniq:
			print ".%s" % c, " ".join([str(i) for i in sol])
	
	print >> sys.stderr, len(solUniq), "solutions"
