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
import utils.myDiags


#############
# FONCTIONS #
#############

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
	print >> sys.stderr, "OK" #, " ".join(["%s:%d" % (anc,len(diagEntry[anc])) for anc in diagEntry])
	


def recCombin(node):

	def getScaffTranslation(anc, i, new):
		if new == None:
			return []
		elif new in listEspeces:
			if i in scaff[new][anc]:
				return scaff[new][anc][i]
			else:
				return []
		else:
			if i in translate[anc][new]:
				j = translate[anc][new][i]
				return [scaff[new][j]]
			else:
				return []
		

	# On arrete la recursion lorsqu'on est sur une espece
	if node in listEspeces:
		return

	# On peut lancer la recursion
	((fils1,_), (fils2,_)) = phylTree.items[node]
	out = outgroup[node]
	print >> sys.stderr, "Recursion sur", node, fils1, fils2, out
	recCombin(fils1)
	recCombin(fils2)
	
	# Initialisation
	lst = diagEntry[node]
	combin = utils.myTools.myCombinator([[i] for i in xrange(len(lst))])
	
	print >> sys.stderr, "Assemblage de", node
	
	# Liste des combinaisons
	for (i1,i2) in utils.myTools.myMatrixIterator(len(lst), len(lst), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	
		#print >> sys.stderr, "diagonale 1", lst[i1]
		#print >> sys.stderr, "diagonale 2", lst[i2]
		
		d1f1 = getScaffTranslation(node, i1, fils1)
		d1f2 = getScaffTranslation(node, i1, fils2)
		d1out = getScaffTranslation(node, i1, out)
		d2f1 = getScaffTranslation(node, i2, fils1)
		d2f2 = getScaffTranslation(node, i2, fils2)
		d2out = getScaffTranslation(node, i2, out)
		
		fils1OK = len(set(d1f1).intersection(d2f1)) > 0
		fils2OK = len(set(d1f2).intersection(d2f2)) > 0
		outOK = len(set(d1out).intersection(d2out)) > 0
		
		fils1OK = fils1OK and (len(set(d1f1).symmetric_difference(d2f1)) == 0)
		fils2OK = fils2OK and (len(set(d1f2).symmetric_difference(d2f2)) == 0)
		outOK = outOK and (len(set(d1out).symmetric_difference(d2out)) == 0)
		
		if (fils1OK and fils2OK) or ((fils1OK or fils2OK) and outOK and options["useOutgroup"]):
		#if (fils1OK and fils2OK):
			combin.addLink([i1,i2])
	
	print >> sys.stderr, "Ecriture de", node
	
	# On rajoute les combinaisons 
	for g in combin:
		if len(g) == 1:
			continue
		nbScaff[node] += 1
		for i in g:
			scaff[node][i] = nbScaff[node]
	
	print >> sys.stderr, "Fin de", node

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("useOutgroup",bool,False), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)

# 2. Chargement et definition des outgroup
genesAnc = {}
diagEntry = {}
outgroup = {}
for anc in phylTree.items:
	diagEntry[anc] = []
	genesAnc[anc] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc, False, False)
	if anc in phylTree.parent:
		freres = phylTree.items[phylTree.parent[anc]]
		for (e,_) in freres:
			if e != anc:
				outgroup[anc] = e
				break
	else:
		outgroup[anc] = None

loadDiagsFile(noms_fichiers["diagsList"], diagEntry)

# Les scaffolds initiaux (au debut: 1 diag = 1 scaff)
# Au fur et a mesure, on fera pointer plusieurs diagonales sur le meme scaffold
print >> sys.stderr, "Initialisation des scaffold ancestraux ..."
scaff = {}
scaff["Opossum"] = {}
nbScaff = {}
for anc in phylTree.items:
	scaff[anc] = [i for i in xrange(len(diagEntry[anc]))]
	nbScaff[anc] = len(scaff[anc])
	scaff["Opossum"][anc] = {}

# La meme structure permet d'associer les diagonales aux chromosomes sur les especes
for esp in listEspeces:
	scaff[esp] = {}
	for anc in phylTree.items:
		scaff[esp][anc] = {}

# On remplit les appartenances des diagonales aux especes
# On cree la table d'association gene-ancestral/diagonales qui le contiennent
print >> sys.stderr, "Initialisation des 'scaffold' modernes et de la table d'association gene/diagonale ..."
positions = {}
for anc in phylTree.items:
	positions[anc] = [[] for i in xrange(len(genesAnc[anc].lstGenes[utils.myGenomes.AncestralGenome.defaultChr]))]
	for i in xrange(len(diagEntry[anc])):
		(_,d,esp) = diagEntry[anc][i]
		for g in d:
			positions[anc][g].append(i)
		for (e,c) in esp:
			scaff[e][anc][i] = scaff[e][anc].get(i, []) + [c]

# Lit toutes les diagonales de anc1 et les mappe sur anc2
def mkTranslation(anc1, anc2):
	sys.stderr.write(".")
	table = {}
	for i in xrange(len(diagEntry[anc1])):
		(_,d,_) = diagEntry[anc1][i]
		# La liste des numeros des diagonales de anc2
		trans = []
		for j in d:
			# La diagonale de anc2 qui correspond a ce gene
			tmp = set([])
			for s in genesAnc[anc1].lstGenes[utils.myGenomes.AncestralGenome.defaultChr][j].names:
				if s in genesAnc[anc2].dicGenes:
					tmp.add(genesAnc[anc2].dicGenes[s][1])
			for s in tmp:
				trans.extend(positions[anc2][s])
		
		# La reponse est la diagonale la plus presente dans trans
		
		# Si elle n'existe pas chez anc2, on arrete
		if len(trans) == 0:
			continue
			
		strans = set(trans)
		ctrans = [(trans.count(x),x) for x in set(trans)]
		ctrans.sort()
		table[i] = ctrans[-1][1]

	return table

# Lit toutes les diagonales de anc1 et les mappe sur anc2
def mkTranslationOutgroup(anc1, anc2):
	sys.stderr.write(".")
	table = {}
	
	#anc = phylTree.getFirstParent(anc1, anc2)
	anc = phylTree.parent[anc1]
	for i in xrange(len(diagEntry[anc1])):
		(_,d,_) = diagEntry[anc1][i]
		# La liste des numeros des diagonales de anc2
		trans = []
		for j in d:
			# La diagonale de anc2 qui correspond a ce gene
			tmp = set([])
			for s in genesAnc[anc1].lstGenes[utils.myGenomes.AncestralGenome.defaultChr][j].names:
				if s in genesAnc[anc].dicGenes:
					tmp.add(genesAnc[anc].dicGenes[s][1])
			tmp2 = set([])
			for k in tmp:
					for s in genesAnc[anc].lstGenes[utils.myGenomes.AncestralGenome.defaultChr][k].names:
						if s in genesAnc[anc2].dicGenes:
							tmp2.add(genesAnc[anc2].dicGenes[s][1])
			for s in tmp2:
				trans.extend(positions[anc2][s])
		
		# La reponse est la diagonale la plus presente dans trans
		
		# Si elle n'existe pas chez anc2, on arrete
		if len(trans) == 0:
			continue
			
		strans = set(trans)
		ctrans = [(trans.count(x),x) for x in set(trans)]
		ctrans.sort()
		table[i] = ctrans[-1][1]

	return table

# La table qui permet de passer d'une diagonale d'un ancetre a celle d'un autre
print >> sys.stderr, "Traduction des diagonales ",
translate = {}
for anc in phylTree.items:
	((fils1,_), (fils2,_)) = phylTree.items[anc]
	translate[anc] = {}
	if fils1 in phylTree.items:
		translate[anc][fils1] = mkTranslation(anc, fils1)
	if fils2 in phylTree.items:
		translate[anc][fils2] = mkTranslation(anc, fils2)
	if outgroup[anc] in phylTree.items:
		translate[anc][outgroup[anc]] = mkTranslationOutgroup(anc, outgroup[anc])
print >> sys.stderr, " OK"
del genesAnc

# 2. Initialisation des outgroup

# On va alterner les merge de 2xFils et les merge de Fils+Outgroup jusqu'a ce
# que plus rien ne change
print >> sys.stderr, "Ini:", [nbScaff[anc] for anc in phylTree.items]
utile = True
while utile:
#for i in xrange(1):
	print >> sys.stderr, "Recursion !"
	length = repr([len(set(scaff[anc])) for anc in phylTree.items])
	recCombin(phylTree.root)
	newLength = repr([len(set(scaff[anc])) for anc in phylTree.items])
	print >> sys.stderr, "OK", newLength
	utile = (length != newLength)


for anc in phylTree.items:
	chrom = set(scaff[anc])
	for c in chrom:
		contenu = set(utils.myMaths.flatten([diagEntry[anc][i][1] for i in xrange(len(diagEntry[anc])) if scaff[anc][i] == c]))
		print "%s\t%d\t%s" % (anc, len(contenu), " ".join([str(x) for x in contenu]))
	

