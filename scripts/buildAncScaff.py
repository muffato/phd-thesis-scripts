#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths


#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, diagEntry):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	for l in f:

		ct = l.split('\t')
		anc = ct[0]
		l = int(ct[1])
		d = [int(x) for x in ct[2].split()]
		esp = [tuple(x.split('/')) for x in ct[3].split()]
		diagEntry[anc].append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

# Lit toutes les diagonales de anc1 et les mappe sur anc2
def mkTranslation(anc1, anc2):
	global genesAnc, diagPositions
	sys.stderr.write(".")
	table = {}
	
	anc = phylTree.parent[anc2]
	for i in xrange(len(diagEntry[anc1])):
		(_,d,_) = diagEntry[anc1][i]
		# La liste des numeros des diagonales de anc2
		trans = []
		for j in d:
			# La diagonale de anc2 qui correspond a ce gene
			tmp = genesAnc[anc].getPosition(genesAnc[anc1].lstGenes[utils.myGenomes.Genome.defaultChr][j])
			tmp2 = set([])
			for (_,k) in tmp:
				tmp2.update(genesAnc[anc2].getPosition(genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr][k]))
			for (_,s) in tmp2:
				trans.extend(diagPositions[anc2][s])
		
		# La reponse est la diagonale la plus presente dans trans
		table[i] = utils.myMaths.unique(trans)
		continue
		
		# Si elle n'existe pas chez anc2, on arrete
		if len(trans) == 0:
			continue
			
		ctrans = [(trans.count(x),x) for x in set(trans)]
		ctrans.sort()
		table[i] = ctrans[-1][1]

	return table


def init():
	global espPositions, genesAnc, diagPositions, translate, combin
	
	# Les scaffolds initiaux (au debut: 1 diag = 1 scaff)
	# Au fur et a mesure, on fera pointer plusieurs diagonales sur le meme scaffold
	print >> sys.stderr, "Initialisation des scaffold ancestraux ..."

	# La meme structure permet d'associer les diagonales aux chromosomes sur les especes
	for esp in listEspeces:
		espPositions[esp] = {}
		for anc in phylTree.items:
			espPositions[esp][anc] = {}

	# On remplit les appartenances des diagonales aux especes
	# On cree la table d'association gene-ancestral/diagonales qui le contiennent
	#print >> sys.stderr, "Initialisation des 'scaffold' modernes et de la table d'association gene/diagonale ..."
	for anc in phylTree.items:
		combin[anc] = utils.myTools.myCombinator([[i] for i in xrange(len(diagEntry[anc]))])
		diagPositions[anc] = [[] for i in xrange(len(genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr]))]
		for i in xrange(len(diagEntry[anc])):
			(_,d,esp) = diagEntry[anc][i]
			for g in d:
				diagPositions[anc][g].append(i)
			for (e,c) in esp:
				if i not in espPositions[e][anc]:
					espPositions[e][anc][i] = set([c])
				else:
					espPositions[e][anc][i].add(c)


	# La table qui permet de passer d'une diagonale d'un ancetre a celle d'un autre
	print >> sys.stderr, "Traduction des diagonales ",
	for anc in phylTree.items:
		((fils1,_), (fils2,_)) = phylTree.items[anc]
		translate[anc] = {}
		if fils1 in phylTree.items:
			translate[anc][fils1] = mkTranslation(anc, fils1)
		if fils2 in phylTree.items:
			translate[anc][fils2] = mkTranslation(anc, fils2)
		if outgroup[anc] in phylTree.items:
			translate[anc][outgroup[anc]] = mkTranslation(anc, outgroup[anc])
	print >> sys.stderr, " OK"


# Fonction qui calcule recursivement les combinaisons sur un noeud de l'arbre
def recCombin(node):

	global espPositions, genesAnc, diagPositions, translate

	vide = set([])
	
	# Retrouve sur quel scaffold de l'ancetre 'new' se trouve la diagonale i de l'ancetre 'anc'
	def getScaffTranslation(anc, i, new):
		if new == None:
			return vide
		elif new in listEspeces:
			if i in espPositions[new][anc]:
				return espPositions[new][anc][i]
			else:
				return vide
		else:
			if i in translate[anc][new]:
				#j = translate[anc][new][i]
				return set([combin[new].dic[j] for j in translate[anc][new][i]])
			else:
				return vide
		

	# On arrete la recursion lorsqu'on est sur une espece
	if node in listEspeces:
		return

	# On peut lancer la recursion
	((fils1,_), (fils2,_)) = phylTree.items[node]
	out = outgroup[node]
	nb = len(diagEntry[node])
	utile = (recCombin(fils1))
	utile = (recCombin(fils2) or utile)
	
	print >> sys.stderr, "Assemblage de %s = [ %s | %s ] ^ %s (%d -> ..." % (node, fils1, fils2, out, combin[node].getNbGrp()),
	
	# Liste des combinaisons
	for (i1,i2) in utils.myTools.myMatrixIterator(nb, nb, utils.myTools.myMatrixIterator.StrictUpperMatrix):
	
		if combin[node].dic[i1] == combin[node].dic[i2]:
			continue

		d1f1 = getScaffTranslation(node, i1, fils1)
		d1f2 = getScaffTranslation(node, i1, fils2)
		d1out = getScaffTranslation(node, i1, out)
		d2f1 = getScaffTranslation(node, i2, fils1)
		d2f2 = getScaffTranslation(node, i2, fils2)
		d2out = getScaffTranslation(node, i2, out)
		
		# Sous-ensembles egaux
		#fils1OK = len(d1f1.intersection(d2f1)) > 0
		#fils2OK = len(d1f2.intersection(d2f2)) > 0
		#outOK = len(d1out.intersection(d2out)) > 0
		
		# Ensembles egaux
		#fils1OK = (d1f1 == d2f1)
		#fils2OK = (d1f2 == d2f2)
		#outOK = (d1out == d2out)
	
		# Ensembles egaux et non vides
		fils1OK = (len(d1f1) > 0) and (d1f1 == d2f1)
		fils2OK = (len(d1f2) > 0) and (d1f2 == d2f2)
		outOK = (len(d1out) > 0) and (d1out == d2out)
	
		if (fils1OK and fils2OK) or ((fils1OK or fils2OK) and outOK and options["useOutgroup"]):
			combin[node].addLink([i1,i2])
			utile = True
	
	print >> sys.stderr, "%d) OK" % (combin[node].getNbGrp())
	return utile


def mergeDiags():
	global combin, diagEntry
	
	for anc in phylTree.items:
		entry = []
		for g in combin[anc]:
			diag = utils.myMaths.unique(utils.myMaths.flatten([diagEntry[anc][i][1] for i in g]))
			esp = utils.myMaths.unique(utils.myMaths.flatten([diagEntry[anc][i][2] for i in g]))
			entry.append( (len(diag),diag,esp) )
		diagEntry[anc] = entry



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("useOutgroup",bool,False), ("mergeDiags",bool,False), ("onlyOnce",bool,False), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root) + ['Opossum']

# Les genes ancestraux
genesAnc = {}
diagEntry = {}
outgroup = {}
espPositions = {}
translate = {}
diagPositions = {}
combin = {}
for anc in phylTree.items:
	diagEntry[anc] = []
	genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc)
	(e1,_) = phylTree.items[anc][0]
	(e2,_) = phylTree.items[anc][1]
	outgroup[e1] = e2
	outgroup[e2] = e1
outgroup[phylTree.root] = None

loadDiagsFile(noms_fichiers["diagsList"], diagEntry)

init()

while recCombin(phylTree.root):

	if options["onlyOnce"]:
		break

	if options["mergeDiags"]:
		# On construit les nouveaux scaffs
		mergeDiags()
		init()


# Impression du resultat
mergeDiags()
for anc in phylTree.items:
	for (l,d,esp) in diagEntry[anc]:
		print "%s\t%d\t%s" % (anc, len(d), " ".join([str(x) for x in d]))

