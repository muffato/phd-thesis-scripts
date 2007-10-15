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


#############
# FONCTIONS #
#############

#
# Calcule la distance inter-genes moyenne entre les deux genes
#   Utilise la fonction d'inference de valeur de phylTree
#
def distInterGenes(tg1, tg2, seuil):

	# Les distances chez chaque espece
	distEsp = {}
	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				# Au dela d'un certain seuil, on considere que l'information n'est plus valable
				if (seuil <= 0) or (x <= seuil):
					# On garde la plus petite distance trouvee
					distEsp[e1] = min(distEsp.get(e1,x), x)
	
	# On fait la liste des especes qui presentent une distance de 1
	lst1Esp = [e for e in distEsp if distEsp[e] == 1]
	
	# On met les 1 dans les noeuds de l'arbre entre les especes
	lst1Anc = set()
	for (e1,e2) in itStrictUpperMatrix(lst1Esp):
		lst1Anc.update(phylTree.dicLinks[e1][e2])
		if anc in lst1Anc:
			return 1
	for a in lst1Anc:
		distEsp[a] = 1

	# En mode outgroups/2 les outgroups qui montrent une distance plus grande que les fils sont supprimes
	if options["useOutgroups"] == 2:
		m = [distEsp[e] for e in ancSpecies if e in distEsp]
		if len(m) > 0:
			m = max(m)
			for e in ancOutgroupSpecies:
				if distEsp.get(e,0) > m:
					del distEsp[e]

	# On calcule par une moyenne les autres distances
	return phylTree.calcDist(distEsp)


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
	if options["useOutgroups"] > 0:
		anc = options["ancestr"]
		while anc in phylTree.parent:
			anc = phylTree.parent[anc]
			# Le genome de l'ancetre superieur
			tmpGenesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
			del tmpGenesAnc.dicGenes
			# Chaque gene
			for g in tmpGenesAnc:
				# Les autres sont les orthologues dans les autres especes
				newGenes = [phylTree.dicGenes[s] for s in g.names if (s in phylTree.dicGenes) and not (s in genesAnc.dicGenes)]
				# On enregistre le lien entre les genes du genome a ordonner et les genes des outgroups
				for (_,i) in genesAnc.getPosition(g.names):
					dicOutgroupGenes[i].update(newGenes)
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
				tmp.extend(dicOutgroupGenes.get( i, []))
				d2.append(tmp)
			genome[c].append(d2)
	del phylTree.dicGenes
	return genome



# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestralDiags", "phylTree.conf"], \
	[("seuilMaxDistInterGenes",float,0), ("nbDecimales",int,2), ("infiniteDist",int,1000000), ("notConstraintPenalty",float,10000), \
	("ancestr",str,""), ("nbConcorde",int,1), ("useOutgroups",int,[0,1,2]), ("withConcordeOutput",bool,False),\
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
diags = loadDiagsFile(noms_fichiers["genomeAncestralDiags"])
genesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[options["ancestr"]])
nbConcorde = max(1, options["nbConcorde"])
mult = pow(10, options["nbDecimales"])
seuil = options["seuilMaxDistInterGenes"]
pen = str(int(mult*options["infiniteDist"]))
add = options["notConstraintPenalty"]
anc = options["ancestr"]
ancSpecies = phylTree.species[anc]
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]

phylTree.initCalcDist(anc, options["useOutgroups"] != 0)
newGenome = rewriteGenome()
concordeInstance = utils.concorde.ConcordeLauncher()

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

		g1 = tab[d1][2*d1-i1]
		g2 = tab[d2][2*d2-i2]
		y = distInterGenes(g1, g2, seuil)
		if y == None:
			return pen
		else:
			return int(mult*y+add)
	
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
		elif i1 > i2:
			# La diagonale est a inverser
			diag.reverse()
			strand.reverse()
			for i in xrange(len(strand)):
				strand[i] *= -1
		# On affiche les genes
		for (i,g) in enumerate(diag):
			print c, strand[i], " ".join(genesAnc.lstGenes[None][g].names)
		
	for sol in lstTot:
		print "# .%s" % c, " ".join([str(i) for i in sol])
	
	print >> sys.stderr, len(lstTot), "solutions OK"

