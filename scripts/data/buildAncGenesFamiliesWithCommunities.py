#! /users/ldog/muffato/python -OO

__doc__ = """
Lit l'arbre phylogenetique et construit les fichiers de genes ancestraux
 pour tous les noeuds en utilisant les fichiers d'orthologues de Biomart.
Affine la liste en exhibant les genes specifiques d'une lignee et ceux qui
 n'ont plus qu'un seul gene moderne correspondant
Utilise le clustering en communautes pour rendre utilisable les releases douteuses
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
import utils.walktrap.myCommunities



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""), ("minRelevance",float,100), ("recursiveConstruction",bool,True), ("homologyFilter",str,"ortholog_one2one,ortholog_one2many,ortholog_many2many"), ("withSubLinks",bool,True), \
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
		for (i1,i2) in utils.myTools.myMatrixIterator(ensGenes, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
			if (i2 in aretes.get(i1,[])) or (i1 in aretes.get(i2,[])):
				nbAretes += 1
		nbAretesAttendu = (len(ensGenes)*(len(ensGenes)-1))/2
		for e in score:
			nbAretesAttendu -= (len(score[e])*(len(score[e])-1))/2
		
		return (nbAretes, nbAretesAttendu, poidsBranches, score)




	# Charge un fichier de definition d'homologues
	def doLoad(s, age):

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
				comb.addLink([gA, gB])

			# Le tableau des aretes entre genes
			if gA not in aretes:
				aretes[gA] = set()
			aretes[gA].add(gB)

		f.close()
		


	# 1. On charge tous les fichiers d'orthologues
	#    On fait les familles par transitivite en ne tenant compte que des liens entre les deux branches
	aretes = {}
	comb = utils.myTools.myCombinator([])
	n = len(phylTree.species[anc])
	print >> sys.stderr, "Construction des familles d'orthologues de %s " % anc,
	for (e1,e2) in utils.myTools.myMatrixIterator(phylTree.species[anc], None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		f = options["orthosFile"] % (phylTree.fileName[e1],phylTree.fileName[e2])
		if phylTree.dicParents[e1][e2] == anc:
			doLoad(f, phylTree.ages[anc])
		elif options["withSubLinks"]:
			doLoad(f, phylTree.ages[anc])
		utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"

	# On cherche a clusteriser les familles transitives
	print >> sys.stderr, "Insertion des genes paralogues intra-especes ",
	for e in phylTree.species[anc]:
		doLoad(options["paras1File"] % phylTree.fileName[e], phylTree.ages[anc]-1)
		utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"
	
	print >> sys.stderr, "Insertion des genes paralogues inter-especes ",
	for (e1,e2) in utils.myTools.myMatrixIterator(phylTree.species[anc], None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		doLoad(options["paras2File"] % (phylTree.fileName[e1],phylTree.fileName[e2]), phylTree.ages[anc])
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

		# A.2/ Tous les liens sont presents
		if nbAretes >= nbAretesAttendu:
			print >> sys.stderr, "Graphe complet"
			clusters = [x]
		
		# A.3/ Dans une des sous-branches, tous les genes sont en 1 copie
		elif min(poidsBranches) == 1:
			print >> sys.stderr, "Branche sans duplication"
			clusters = [x]

		# B. On est oblige de clusteriser
		else:
			# B.1/ On lance les communautes
			lstCommunitiesOrig = utils.walktrap.myCommunities.launchCommunitiesBuild(items = x, edgesDict = aretes)[0]
			
			# B.2/ On selectionne celles qui sont convenables
			lstCommunities = []
			for comm in lstCommunitiesOrig:
			
				print >> sys.stderr, "Test de [alpha=%f relevance=%f parts=%d N/A=%d/%d] :" % (comm[0],comm[1],len(comm[2]),len(comm[3]),len(x)),
				(alpha,relevance,clusters,lonely) = comm
				if options["graphDirectory"] != "" and len(lonely) == 0:
					fa = utils.myTools.myOpenFile(options["graphDirectory"] + '/graph-%f-%f-%d-%d' % (relevance,alpha,len(clusters),len(lonely)), 'w')
					print >> fa, "graph {"
					interv = 85
					for ci in xrange(1,len(clusters)+1):
						r = ci % 4
						g = (ci / 4) % 4
						b = (ci / 16) % 4
						for cc in clusters[ci-1]:
							print >> fa, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (cc, 85*r,85*g,85*b)
					for (i1,i2) in utils.myTools.myMatrixIterator(x, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
						if i2 in aretes.get(i1, []):
							ss = 1
						elif i1 in aretes.get(i2, []):
							ss = 1
						elif phylTree.dicGenes[i1][0] == phylTree.dicGenes[i2][0]:
							ss = -1
						else:
							ss = 0
						if ss == -1:
							print >> fa, "%s -- %s [style=\"dotted\"]" % (i1, i2)
						elif ss == 0.5:
							print >> fa, "%s -- %s [style=\"dashed\"]" % (i1, i2)
						elif ss == 1:
							print >> fa, "%s -- %s" % (i1, i2)
					print >> fa, "}"
					fa.close()


				# Ne nous interessent que les clusterisations completes
				if len(comm[-1]) != 0:
					print >> sys.stderr, "genes oublies"
					continue
				
				# Une clusterisation est correcte si tous les clusters sont repartis sur les deux sous-branches
				tmpNbEsp = []
				for c in comm[2]:
					(nbAr, nbArAtt, poidsBranches, score) = getStats(c)
					tmpNbEsp.append(   float(len([e for e in score if len(score[e]) > 0])) )
					if len([e for e in poidsBranches if e > 0]) == 1:
						print >> sys.stderr, "Clusterisation incoherente"
						break
				else:
					print >> sys.stderr, min(tmpNbEsp),
					# Ne nous interessent que les clusterisations avec une relevance suffisante
					if comm[1] >= options["minRelevance"]:
						print >> sys.stderr, "Communaute recevable"
						lstCommunities.append(comm)
					else:
						print >> sys.stderr, "relevance insuffisante"
					
			# On trie suivant la meilleure relevance
			lstCommunities.sort(key = operator.itemgetter(1), reverse = True)

			# B.3/ Aucune clusterisation satisfaisante
			if len(lstCommunities) == 0:
				print >> sys.stderr, "Aucune communaute satisfaisante"
				clusters = [x]
				
			# B.4/ Unique solution
			elif len(lstCommunities) == 1:
				clusters = lstCommunities[0][2]
				print >> sys.stderr, "Resultat unique [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
			# B.5/ Le choix s'impose
			else:
				clusters = lstCommunities[0][2]
				print >> sys.stderr, "Hesitation ... Choix: [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
			
		# Ecriture
		for c in clusters:
			res.addLink(c)
	
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

