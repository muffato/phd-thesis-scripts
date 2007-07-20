#! /users/ldog/muffato/python -OO

__doc__ = """
Lit l'arbre phylogenetique et construit les fichiers de genes ancestraux
 pour tous les noeuds en utilisant les fichiers d'orthologues de Biomart.
Affine la liste en exhibant les genes specifiques d'une lignee et ceux qui
 n'ont plus qu'un seul gene moderne correspondant
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import operator
import utils.myGenomes
import utils.myTools
import utils.myPhylTree



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""), ("recursiveConstruction",bool,True), ("homologyFilter",str,"ortholog_one2one,ortholog_one2many,ortholog_many2many"), ("withSubLinks",bool,True), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("one2oneFile",str,"~/work/data/ancGenes/one2one.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("paras2File",str,"~/work/data/orthologs/paras.%s.%s.list.bz2"), \
	("paras1File",str,"~/work/data/paralogs/paras.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince(options["ancestr"], options["genesFile"])
homologies = set(options["homologyFilter"].split(","))

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
			if champs[-1] in homologies:
				combin.addLink([gA, gB])

		f.close()
		


	# 1. On charge tous les fichiers d'orthologues
	#    On fait les familles par transitivite en ne tenant compte que des liens entre les deux branches
	aretes = {}
	comb = utils.myTools.myCombinator([])
	n = len(phylTree.species[anc])
	print >> sys.stderr, "Construction des familles d'orthologues de %s " % anc,
	for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.species[anc]):
		f = options["orthosFile"] % (phylTree.fileName[e1],phylTree.fileName[e2])
		if phylTree.dicParents[e1][e2] == anc:
			doLoad(f, phylTree.ages[anc], comb)
		elif options["withSubLinks"]:
			doLoad(f, phylTree.ages[anc], comb)
		utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"

	# On cherche a clusteriser les familles transitives
	if "within_species_paralog" in homologies:
		print >> sys.stderr, "Insertion des genes paralogues intra-especes ",
		for e in phylTree.species[anc]:
			doLoad(options["paras1File"] % phylTree.fileName[e], phylTree.ages[anc]-1, comb)
			utils.myTools.stderr.write('.')
		print >> sys.stderr, " OK"
		
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "..."
	
	res = utils.myTools.myCombinator([])
	for x in comb:
	
		(nbAretes,nbAretesAttendu,poidsBranches,score) = getStats(x)

		# A.1/ On elimine les familles specifiques d'une sous-branche
		if len([e for e in poidsBranches if e > 0]) == 1:
			#print >> sys.stderr, "Sous-famille"
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
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % phylTree.fileName[anc], 'w')
	ff = utils.myTools.myOpenFile(options["one2oneFile"] % phylTree.fileName[anc], 'w')
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

