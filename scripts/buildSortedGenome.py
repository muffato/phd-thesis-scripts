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

matIt = utils.myTools.myMatrixIterator
matItmWhole = utils.myTools.myMatrixIterator.WholeMatrix
matItmSUpper = utils.myTools.myMatrixIterator.StrictUpperMatrix


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
	for ((e1,c1,i1), (e2,c2,i2)) in matIt(tg1, tg2, matItmWhole):
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
	for (e1,e2) in matIt(lst1Esp, None, matItmSUpper):
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


#
# Reecrit le genome a trier comme suite des positions des genes dans les genomes modernes
#
def rewriteGenome():

	# On etend la liste des genes ancestraux pour utiliser les outgroup en remontant l'arbre jusqu'a la racine
	dicOutgroupGenes = {}
	anc = options["ancestr"]
	while anc in phylTree.parent:
		anc = phylTree.parent[anc]
		# Le genome de l'ancetre superieur
		tmpGenesAnc = utils.myGenomes.loadGenome(options["ancGenesFile"] % phylTree.fileName[anc])
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
	# Chaque gene devient la liste des distances inter-genes chez chaque espece
	genome = {}
	for c in genesAnc.lstChr:
		genome[c] = []
		for (i,g) in enumerate(genesAnc.lstGenes[c]):
			tmp = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
			tmp.extend(dicOutgroupGenes.get( (c,i), []))
			genome[c].append(tmp)
	del phylTree.dicGenes
	return genome


#
# Ecrit la matrice du chromosome c et lance concorde
#
def sortChromosome(c):

	tab = genome[c]
	n = len(tab)

	mult = pow(10, options["nbDecimales"])
	seuil = options["seuilMaxDistInterGenes"]
	pen = str(int(mult*options["penalite"]))
	
	print >> sys.stderr
	print >> sys.stderr, "- Chromosome %s (%d genes) -" % (c, n)
	
	# La matrice des distances intergenes
	print >> sys.stderr, "Ecriture de la matrice ... ",
	f = open(nom, 'w', 65536)
	print >> f, "NAME: CHRANC"
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	print >> f, "0 " * n
	for (i,t) in enumerate(tab):
		for j in xrange(i+1,n):
			y = distInterGenes(t, tab[j], seuil)
			if y == None:
				print >> f, pen,
			elif y == 1:
				print >> f, "0",
			else:
				print >> f, str(int(mult*y)),
	print >> f, "EOF"
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	lstTot = []
	for i in xrange(nbConcorde):
		comm = options["concordeExec"] + ' -x ' + nom
		if options["withConcordeOutput"]:
			os.system(comm + ' >&2')
		else:
			os.system(comm + ' > /dev/null')
		if os.access(nom + ".sol", os.R_OK):
			lstTot.append(ConcordeFile(nom + ".sol"))
		os.system('rm -f 0%s* %s*' % (nom,nom) )
		sys.stderr.write(".")

	# On remet chaque liste dans le meme sens que la premiere
	for t in lstTot[1:]:
		if not t.isMemeSens(lstTot[0]):
			t.reverse()

	if len(lstTot) == 0:
		for g in genesAnc.lstGenes[c]:
			print c,
			if (options["nbConcorde"] > 1) and options["withConcordeStats"]:
				print 1,
			print " ".join(g.names)
	else:
		for (i,g) in enumerate(lstTot[0].res):
			print c,
			if (options["nbConcorde"] > 1) and options["withConcordeStats"]:
				print len(set([s.res[i] for s in lstTot])),
			print " ".join(genesAnc.lstGenes[c][g-1].names)
	
	solUniq = utils.myMaths.unique([l.res for l in lstTot])
	print >> sys.stderr, len(solUniq), "solutions"
	if options["withConcordeStats"]:
		for sol in solUniq:
			print ".%s" % c, " ".join([str(i) for i in sol])


	
#####################################################
# Cette classe gere un fichier resultat de Concorde #
#####################################################
class ConcordeFile:

	def __init__(self, nom):
		tmp = []
		f = utils.myTools.myOpenFile(nom, 'r')
		for ligne in f:
			for x in ligne.split():
				tmp.append(int(x))
		
		f.close()
		i = tmp.index(0)
		self.res = tmp[i+1:] + tmp[1:i]
		self.buildIndex()
			

	def buildIndex(self):
		self.dic = {}
		for (i,x) in enumerate(self.res):
			self.dic[x] = i
	
	def isMemeSens(self, other):
		s = 0
		for i in xrange(len(self.res)-1):
			if other.dic[self.res[i]] < other.dic[self.res[i+1]]:
				s += 1
			else:
				s -= 1
		return s>0

	def reverse(self):
		self.res.reverse()
		self.buildIndex()
	


# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestral", "phylTree.conf"], \
	[("ancestr",str,""), ("seuilMaxDistInterGenes",float,0), ("nbDecimales",int,2), ("penalite",int,1000000), \
	("useOutgroups",int,[0,1,2]), ("nbConcorde",int,-1), ("withConcordeOutput",bool,False), ("withConcordeStats",bool,False),\
	("concordeExec",str,"~/work/scripts/concorde"),\
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
phylTree.loadAllSpeciesSince(None, options["genesFile"])
del phylTree.dicGenomes
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genomeAncestral"])
genome = rewriteGenome()
nom = "mat%08d" % ((os.getpid() ^ os.getppid() ^ random.randint(1,16777215)) & 16777215)
nbConcorde = max(1, options["nbConcorde"])
anc = options["ancestr"]
ancSpecies = phylTree.species[anc]
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]
phylTree.initCalcDist(anc, options["useOutgroups"] != 0)

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:
	sortChromosome(c)

os.system('rm -f *%s*' % nom )

