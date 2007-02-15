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
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myCommunities
import utils.myPsOutput




########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""),("minRelevance",float,0.3), ("recursiveConstruction",bool,True), ("graphDirectory",str,""), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/full/genes.%s.list.bz2"), \
	("one2oneFile",str,"~/work/data/one2one/one2one.%s.list.bz2"), \
	("orthoFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince(options["ancestr"], options["genesFile"])
utils.myPsOutput.initColor()


def buildAncFile(anc, lastComb):

	def doLoad(s):

		f = utils.myTools.myOpenFile(s, 'r')
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			gA = champs[0]
			gB = champs[3]
			
			# Si les genes ont ete separes au noeud superieur, on les laisse separes
			if lastComb.dic.get(gA,0) != lastComb.dic.get(gB,0):
				continue
			
			comb.addLink([gA, gB])
			if gA not in aretes:
				aretes[gA] = dict([])
			
			# Les apparent_one2one comptent pour 1/2 dans les poids des aretes
			if champs[6].startswith("apparent"):
				aretes[gA][gB] = 0.5
				aretes[gA][gB] = 1
			else:
				aretes[gA][gB] = 1

		f.close()
		


	# 1. on combine tous les fichiers d'orthologues
	esp = [[phylTree.commonNames[x][0] for x in s] for s in phylTree.branchesSpecies[anc]]
	aretes = {}
	comb = utils.myTools.myCombinator([])
	
	print >> sys.stderr, "Construction des familles d'orthologues de", anc, ":", "-".join(utils.myMaths.flatten(esp)), "",
	for (i,j) in utils.myTools.myMatrixIterator(len(esp), len(esp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		for e1 in esp[i]:
			for e2 in esp[j]:
				doLoad(options["orthoFile"] % (e1,e2))
				utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"
	
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "..."
	
	s = anc.replace('/', '_').replace(' ', '_')
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % s, 'w')
	ff = utils.myTools.myOpenFile(options["one2oneFile"] % s, 'w')
	nbA = 0
	nbO = 0
	res = utils.myTools.myCombinator([])
	for x in comb:
	
		score = phylTree.findFamilyComposition(x)
		poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]

		nbAretes = 0
		for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			if x[i1] in aretes and x[i2] in aretes[x[i1]]:
				nbAretes += 1
			elif x[i2] in aretes and x[i1] in aretes[x[i2]]:
				nbAretes += 1
		nbAretesAttendu = 0
		tmp = [sum([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
		for (i,j) in utils.myTools.myMatrixIterator(len(tmp), len(tmp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			nbAretesAttendu += tmp[i]*tmp[j]

		# A. Des conditions qui aident
		# A.1/ Dans une des sous-branches, tous les genes sont en 1 copie
		if nbAretes == nbAretesAttendu:
			print >> sys.stderr, "Graphe complet"
			clusters = [x]
		
		# A.2/ Dans une des sous-branches, tous les genes sont en 1 copie
		elif min(poidsBranches) == 1:
			print >> sys.stderr, "Branche sans duplication"
			clusters = [x]

		# B. On est oblige de clusteriser
		else:
		
			# B.1/ On lance les communautes
			lstCommunitiesOrig = utils.myCommunities.launchCommunitiesBuild1(x, aretes)
			
			# B.2/ On selectionne celles qui sont convenables
			lstCommunities = []
			for comm in lstCommunitiesOrig:
			
				print >> sys.stderr, "Test de [alpha=%f relevance=%f parts=%d N/A=%d/%d] :" % (comm[0],comm[1],len(comm[2]),len(comm[3]),len(x)),
				(alpha,relevance,clusters,lonely) = comm
				if len(options["graphDirectory"]) > 0:
					fa = open(options["graphDirectory"] + '/graph-%f-%f-%d-%d' % (relevance,alpha,len(clusters),len(lonely)), 'w')
					print >> fa, "graph {"
					for ci in xrange(len(clusters)):
						c = clusters[ci]
						(r,g,b) = utils.myPsOutput.colorTable[utils.myPsOutput.color[str(ci+1)]]
						for cc in c:
							print >> fa, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (cc, int(255*r),int(255*g),int(255*b))
					for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
						if x[i1] in aretes and x[i2] in aretes[x[i1]]:
							ss = aretes[x[i1]][x[i2]]
						elif x[i2] in aretes and x[i1] in aretes[x[i2]]:
							ss = aretes[x[i2]][x[i1]]
						elif phylTree.dicGenes[x[i1]][0] == phylTree.dicGenes[x[i2]][0]:
							ss = -1
						else:
							ss = 0
						if ss == -1:
							print >> fa, "%s -- %s [style=\"dotted\"]" % (x[i1], x[i2])
						elif ss == 0.5:
							print >> fa, "%s -- %s [style=\"dashed\"]" % (x[i1], x[i2])
						elif ss == 1:
							print >> fa, "%s -- %s" % (x[i1], x[i2])
					print >> fa, "}"
					fa.close()


				# Ne nous interessent que les clusterisations completes
				if len(comm[-1]) != 0:
					print >> sys.stderr, "genes oublies"
					continue
				
				# Une clusterisation est correcte si tous les clusters sont repartis sur les deux sous-branches
				for c in comm[2]:
					score = phylTree.findFamilyComposition(c)
					poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
					if len([e for e in poidsBranches if e > 0]) == 1:
						print >> sys.stderr, "Clusterisation incoherente"
						break
				else:
			
					# Ne nous interessent que les clusterisations avec une relevance suffisante
					if comm[1] >= options["minRelevance"]:
						print >> sys.stderr, "Communaute recevable"
						lstCommunities.append(comm)
					else:
						print >> sys.stderr, "relevance insuffisante"
						continue
					
			# On trie suivant la meilleure relevance
			lstCommunities.sort(key = operator.itemgetter(1), reverse = True)

			# B.3/ Aucune clusterisation satisfaisante
			if len(lstCommunities) == 0:
				print >> sys.stderr, "Aucune communaute satisfaisante"
				clusters = [x]
				
			# B.4/ Unique solution
			elif len(lstCommunities) == 1:

				# TODO: Verifier le score de relevance !!
				clusters = lstCommunities[0][2]
				print >> sys.stderr, "Resultat unique [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
			
			# B.5/ Le choix s'impose
			else:
				clusters = lstCommunities[0][2]
				print >> sys.stderr, "Hesitation ... Choix: [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
			
			
		nbA += len(clusters)

		# Ecriture
		for c in clusters:
			res.addLink(c)
			print >> f, " ".join(c)
			
			score = dict( [(e,[]) for e in phylTree.dicGenomes] )
			for i in c:
				score[phylTree.dicGenes[i][0]].append(i)
			l = [score[e][0] for e in score if len(score[e]) == 1]
			if len(l) >= 1:
				print >> ff, " ".join(l)
				nbO += 1
	
	comb = None

	add = 0
	# 3. On rajoute les familles heritees
	print >> sys.stderr, "Rajout des familles heritees ...",
	fils = set(phylTree.species[anc])
	for gr in lastComb:
		oublies = [g for g in gr if (g not in res.dic) and (phylTree.dicGenes[g][0] in fils)]
		if len(oublies) > 0:
			print >> f, " ".join(oublies)
			nbA += 1
			print >> ff, " ".join(oublies)
			nbO += 1
			add += 1
			res.addLink(oublies)

	f.close()
	ff.close()
	
	print >> sys.stderr, "OK (%d\\%d+%d)" % (nbA,nbO,add)

	if options["recursiveConstruction"]:
		for (esp,_) in phylTree.items[anc]:
			if esp not in phylTree.dicGenomes:
				buildAncFile(esp, res)

buildAncFile(options["ancestr"], utils.myTools.myCombinator([]))
