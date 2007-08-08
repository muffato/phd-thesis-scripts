#! /users/ldog/muffato/python -OO

__doc__ = """
Lit l'arbre phylogenetique et construit les fichiers de genes ancestraux
 pour tous les noeuds en utilisant les fichiers d'orthologues de Biomart.
Affine la liste en exhibant les genes specifiques d'une lignee et ceux qui
 n'ont plus qu'un seul gene moderne correspondant
"""

# Librairies
import os
import sys
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""), ("recursiveConstruction",bool,True), ("OUT.directory",str,""), \
	("OUT.ancGenesFile",str,"ancGenes/ancGenes.%s.list.bz2"), \
	("OUT.one2oneFile",str,"ancGenes/one2one.%s.list.bz2"), \
	("IN.genesFile",str,"genes/genes.%s.list.bz2"), \
	("IN.orthosFile",str,"orthologs/orthos.%s.%s.list.bz2"), \
	], \
	__doc__ \
)
	

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTancGenesFile = os.path.join(options["OUT.directory"], options["OUT.ancGenesFile"])
OUTone2oneFile = os.path.join(options["OUT.directory"], options["OUT.one2oneFile"])
INgenesFile = os.path.join(options["OUT.directory"], options["IN.genesFile"])
INorthosFile = os.path.join(options["OUT.directory"], options["IN.orthosFile"])
for dir in [OUTancGenesFile, OUTone2oneFile, INgenesFile, INorthosFile]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince(options["ancestr"], INgenesFile)

def buildAncFile(anc, lastComb):


	# Analyse une composante connexe: renvoie le nombre d'aretes et le nombre attendu ainsi que le detail du nombre de genes pour chaque espece
	def getStats(ensGenes):
		score = phylTree.findFamilyComposition(ensGenes)
		poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]

		nbAretes = 0
		for (i1,i2) in utils.myTools.myIterator.tupleOnStrictUpperList(ensGenes):
			if (i2 in aretes.get(i1,[])) or (i1 in aretes.get(i2,[])):
				nbAretes += 1
		nbAretesAttendu = (len(ensGenes)*(len(ensGenes)-1))/2
		for e in score:
			nbAretesAttendu -= (len(score[e])*(len(score[e])-1))/2
		
		return (nbAretes, nbAretesAttendu, poidsBranches, score)




	# Charge un fichier de definition d'homologues
	def doLoad(s, age, combin):

		f = utils.myTools.myOpenFile(s, 'r')
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne[:-1].split('\t')
			gA = champs[0]
			gB = champs[3]
			
			# Si les genes ont ete separes au noeud superieur, on les laisse separes
			if (gA in lastComb.dic) and (gB in lastComb.dic) and (lastComb.dic[gA] != lastComb.dic[gB]):
				continue

			# Si la date du lien d'orthologie est trop ancienne, on ne le prend pas en compte
			if phylTree.ages[champs[6]] > age:
				continue

			# Si on ne veut pas les apparent par exemple
			if champs[-1] in ["ortholog_one2one","ortholog_one2many","ortholog_many2many"]:
				combin.addLink([gA, gB])

		f.close()
		


	# 1. On charge tous les fichiers d'orthologues
	#    On fait les familles par transitivite en ne tenant compte que des liens entre les deux branches
	aretes = {}
	comb = utils.myTools.myCombinator([])
	n = len(phylTree.species[anc])
	print >> sys.stderr, "Construction des familles d'orthologues de %s " % anc,
	for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.species[anc]):
		if phylTree.dicParents[e1][e2] == anc:
			f = INorthosFile % (phylTree.fileName[e1],phylTree.fileName[e2])
			doLoad(f, phylTree.ages[anc], comb)
			utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"

	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "..."
	
	res = utils.myTools.myCombinator([])
	for x in comb:
	
		(nbAretes,nbAretesAttendu,poidsBranches,score) = getStats(x)

		# A.1/ On elimine les familles specifiques d'une sous-branche
		if len([e for e in poidsBranches if e > 0]) == 1:
			print >> sys.stderr, "Sous-famille"
			continue

		res.addLink(x)
	
	comb = None
	aretes = None

	add = 0
	# 3. On rajoute les familles heritees
	print >> sys.stderr, "Rajout des familles heritees ...",
	fils = set(phylTree.species[anc])
	for gr in lastComb:
		genesOK = [g for g in gr if phylTree.dicGenes[g][0] in fils]
		deja = set([res.dic[g] for g in genesOK if g in res.dic])
		oublies = [g for g in genesOK if g not in res.dic]
		if len(oublies) == 0:
			continue
		# Si les genes oublies iraient bien dans 1 famille, on les y insere
		if len(deja) == 1:
			res.addLink(genesOK)
		# Sinon, on cree une nouvelle famille
		else:
			res.addLink(oublies)
			add += 1

	# On ecrit les fichiers
	f = utils.myTools.myOpenFile(OUTancGenesFile % phylTree.fileName[anc], 'w')
	ff = utils.myTools.myOpenFile(OUTone2oneFile % phylTree.fileName[anc], 'w')
	nbA = 0
	nbO = 0
	for c in res:
		print >> f, " ".join(c)
		nbA += 1
		score = phylTree.findFamilyComposition(c)
		l = [score[e][0] for e in score if len(score[e]) == 1]
		print >> ff, " ".join(l)
		if len(l) >= 1:
			nbO += 1
	f.close()
	ff.close()
	
	print >> sys.stderr, "OK (%d\\%d+%d)" % (nbA,nbO,add)

	if options["recursiveConstruction"]:
		for esp in phylTree.branches[anc]:
			if esp not in phylTree.dicGenomes:
				buildAncFile(esp, res)

buildAncFile(options["ancestr"], utils.myTools.myCombinator([]))

