#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import os
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree


# FONCTIONS #

def distInterGenes(tg1, tg2):

	# Les distances chez chaque espece
	distEsp = {}
	for ((e1,c1,i1), (e2,c2,i2)) in utils.myTools.myMatrixIterator(tg1, tg2, utils.myTools.myMatrixIterator.WholeMatrix):
		if e1 == e2 and c1 == c2:
			x = abs(i1-i2)
			# Au dela d'un certain seuil, on considere que l'information n'est plus valable
			if (options["seuilMaxDistInterGenes"] <= 0) or (x <= options["seuilMaxDistInterGenes"]):
				# On garde la plus petite distance trouvee
				distEsp[e1] = min(distEsp.get(e1,x), x)
	
	# On fait la liste des especes qui presentent une distance de 1
	lst1Esp = [e for e in distEsp if distEsp[e] == 1]
	
	# On met les 1 dans les noeuds de l'arbre entre les especes
	lst1Anc = set()
	for (e1,e2) in utils.myTools.myMatrixIterator(lst1Esp, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		lst1Anc.update(phylTree.dicLinks[e1][e2])
		if anc in lst1Anc:
			return 1
	for a in lst1Anc:
		distEsp[a] = 1

	# En mode outgroups/2 les outgroups qui montrent une distance plus grande que les fils sont supprimes
	if options["useOutgroups"] == 2:
		m = [distEsp[e] for e in ancSpecies if e in distEsp]
		if len(m) > 0:
			eNO = [e for e in ancOutgroupSpecies if distEsp.get(e,0) > max(m)]
			for e in eNO:
				del distEsp[e]

	# On calcule par une moyenne les autres distances
	return phylTree.calcDist(distEsp)

	
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
		#i = tmp.index(0)
		#self.res = tmp[i+1:] + tmp[1:i]
		self.res = tmp[1:]
	


# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestral", "phylTree.conf"], \
	[("ancestr",str,""), ("seuilMaxDistInterGenes",float,0), ("nbDecimales",int,2), ("penalite",int,1000000), \
	("useOutgroups",int,[0,1,2]), ("withConcordeOutput",bool,False), \
	("concordeExec",str,"~/work/scripts/concorde"),\
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	"Trie les gens dans l'ordre indique par l'arbre phylogenetique" \
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

# On etend la liste des genes ancestraux pour utiliser les outgroup
anc = options["ancestr"]
ancSpecies = phylTree.species[anc]
ancOutgroupSpecies = phylTree.outgroupSpecies[anc]

dicOutgroupGenes = {}
# On remonte l'arbre jusqu'a la racine
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

# On reecrit le genome en terme d'especes
genome = {}
for c in genesAnc.lstChr:
	genome[c] = []
	for i in xrange(len(genesAnc.lstGenes[c])):
		g = genesAnc.lstGenes[c][i]
		tmp = [phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes]
		tmp.extend(dicOutgroupGenes.get( (c,i), []))
		genome[c].append(tmp)
del phylTree.dicGenes
del dicOutgroupGenes

nom = "mat%08d" % ((os.getpid() ^ os.getppid() ^ random.randint(1,16777215)) & 16777215)
mult = pow(10, options["nbDecimales"])
phylTree.initCalcDist(options["ancestr"], options["useOutgroups"] != 0)

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:

	tab = genome[c]
	n = len(tab)
	
	print >> sys.stderr
	print >> sys.stderr, "- Chromosome %s (%d genes) -" % (c, n)
	
	mat = [range(n+1) for i in xrange(n+1)]
	
	print >> sys.stderr, "Ecriture de la matrice ... ",
	f = open(nom, 'w')
	print >> f, "NAME: CHRANC"
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	print >> f, "0 " * n
	for i in xrange(n):
		for j in xrange(i+1,n):
			y = distInterGenes(tab[i], tab[j])
			mat[i+1][j+1] = mat[j+1][i+1] = y
			if y == None:
				print >> f, int(mult*options["penalite"]),
			elif y == 1:
				print >> f, 0,
			else:
				print >> f, int(mult*y),
		print >> f
	for i in xrange(n+1):
		mat[i][0] = mat[0][i] = 0

	print >> f, "EOF"
	f.close()

	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	
	comm = options["concordeExec"] + ' -x ' + nom
	if options["withConcordeOutput"]:
		os.system(comm + ' >&2')
	else:
		os.system(comm + ' > /dev/null')
	if os.access(nom + ".sol", os.R_OK):
		res = ConcordeFile(nom + ".sol").res
	else:
		res = None
		print >> sys.stderr, "No solution found"
	os.system('rm -f 0%s* %s*' % (nom,nom) )
	sys.stderr.write(".")

	if res == None:
		continue

	#for i in xrange(n):
	#	print c," ".join(genesAnc.lstGenes[c][i].names)
	
	reversalPoints = []
	tmp = [(i,res[i],res[i+1],mat[res[i]][res[i+1]]) for i in xrange(n)]
	for (i,xi,xip,diip) in tmp:
		for (j,xj,xjp,djjp) in tmp[:i-1]:
			dij = mat[xi][xj]
			if dij == None:
				continue
			dijpp = mat[xip][xjp]
			if dijpp == None:
				continue
			if int(mult * (diip+djjp-dij-dijpp)) == 0:
				# On a detecte une inversion potentielle
				reversalPoints.append( (i,j) )
	
	print >> sys.stderr, len(reversalPoints), "inversions possibles"


os.system('rm -f *%s*' % nom )


