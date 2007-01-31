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
	["genesList.conf", "phylTree.conf"], \
	[("homologyLevels",str,"ortholog_one2many,ortholog_many2many,apparent_ortholog_one2one,ortholog_one2one"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("one2oneFile",str,"~/work/data/one2one/one2one.%s.list.bz2"), \
	("orthoFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.listSpecies)
homologies = options["homologyLevels"].split(",")
utils.myPsOutput.initColor()


def buildAncFile(anc, lastComb):

	def doLoad(s):

		f = utils.myTools.myOpenFile(s, 'r')
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			gA = champs[0]
			gB = champs[3]
			
			if gA not in aretes:
				aretes[gA] = dict([])
			aretes[gA][gB] = 1
			
			# Les apparent_one2one comptent pour 1/2 dans les poids des aretes
			if champs[6].startswith("apparent"):
				aretes[gA][gB] = 0.5
			comb.addLink([gA, gB])

		f.close()
		

	def hasTrueLinksBetweenClusters(x, (alpha,relevance,clusters,lonely)):
	
		for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			realLink = (((i1,i2) in tmpAretes) and ((i1,i2) not in tmpAretesApparent))
			
			memeCluster = False
			for c in clusters:
				if (x[i1] in c) and (x[i2] in c):
					memeCluster = True
					break
			if realLink and (not memeCluster):
				print >> sys.stderr, "break", x[i1], x[i2]
				return True
		return False
	

	# 1. on combine tous les fichiers d'orthologues
	comb = utils.myTools.myCombinator([])
	esp = phylTree.species[anc]
	aretes = {}
	
	print >> sys.stderr, "Construction des familles d'orthologues de", anc, ":", "-".join(esp), "",
	for (i,j) in utils.myTools.myMatrixIterator(len(esp), len(esp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		doLoad(options["orthoFile"] % (esp[i],esp[j]))
		utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"
	
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "..."
	
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % anc, 'w')
	ff = utils.myTools.myOpenFile(options["one2oneFile"] % anc, 'w')
	nbA = 0
	nbO = 0
	res = set([])
	for x in comb:
	
		# A. Composition du gene ancestral suivant les familles
		score = utils.myGenomes.findFamilyComposition(x)
		poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
		
		# B. On enleve les familles specifiques d'une unique sous-branche et non heritees du noeud superieur
		if (x[0] not in lastComb) and (len([e for e in poidsBranches if e > 0]) == 1):
			continue

		# C. On calcule les aretes du graphe
		tmpAretes = []
		tmpAretesApparent = []
		tmpAretesMemeEspece = []
		for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			g1 = x[i1]
			g2 = x[i2]
			if (g1 in aretes.get(g2,[])) or (g2 in aretes.get(g1,[])):
				try:
					t1 = aretes[g1][g2]
				except Exception:
					t1 = aretes[g2][g1]
				if t1 == 0.5:
					tmpAretesApparent.append( (i1,i2) )
				tmpAretes.append( (i1,i2) )
			if (geneBank.dicGenes[g1][0] == geneBank.dicGenes[g2][0]):
				tmpAretesMemeEspece.append( (i1,i2) )
				#if g1 not in aretes:
				#	aretes[g1] = dict([])
				#aretes[g1][g2] = 0.5

		
		# D. Filtres
		# D.1/ Graphe complet -> 1 cluster
		if (len(tmpAretes) + len(tmpAretesMemeEspece)) == (len(x)*(len(x)-1))/2:
			print >> sys.stderr, "Graphe complet"
			clusters = [x]
			
		# D.2/ Tous les genes en 1 copie dans une des branches -> 1 cluster
		elif min(poidsBranches) == 1:
			print >> sys.stderr, "Branche sans duplication"
			clusters = [x]

		# D.3/ On est oblige de clusteriser
		else:
			#(relev,clusters) = utils.myCommunities.launchCommunitiesBuild(len(x), test, minCoverage=0.9, minRelevance=0.4)
			#(relev,clusters) = utils.myCommunities.launchCommunitiesBuild(len(x), test, bestRelevance = False)
			lstCommunitiesOrig = utils.myCommunities.launchCommunitiesBuild1(x, aretes)
			
			# Ne nous interessent que les clusterisations completes triees dans l'ordre des fusions
			lstCommunities = [c for c in lstCommunitiesOrig if len(c[-1]) == 0]
			lstCommunities.sort()
			lstCommunities.sort(key = operator.itemgetter(1), reverse = True)

			#lstCommunities = [c for c in lstCommunitiesOrig if isO2Ocompliant(x, c)]
			#if len([c for c in lstCommunitiesOrig if (not isO2Ocompliant(x, c)) and (len(c[-1]) == 0) and (c[1] >= 0.4)]) >= 1:
			#	print >> sys.stderr, "souci"
			
			# D.3.a/ Si il reste toujours des genes non assignes -> 1 cluster
			if len(lstCommunities) == 0:
				print >> sys.stderr, "Aucune communaute satisfaisante"
				clusters = [x]
				
			# D.3.b/ Cluster unique sans equivoque
			elif len(lstCommunities) == 1:

				# TODO: Verifier le score de relevance !!
				clusters = lstCommunities[0][2]
				print >> sys.stderr, "Resultat unique [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
			
			# D.3.c/ Le choix s'impose
			else:
				print >> sys.stderr, "Hesitations !"
			
			for comm in lstCommunities:
				if not hasTrueLinksBetweenClusters(x, comm):
					print >> sys.stderr, "Winner", comm
					break
			else:
				print >> sys.stderr, "Nothing"
			
			for (alpha,relevance,clusters,_) in lstCommunities:
				fa = open('/users/ldog/muffato/work/tutu/graph-%f-%f' % (relevance,alpha), 'w')
				print >> fa, "graph {"
				for ci in xrange(len(clusters)):
					c = clusters[ci]
					(r,g,b) = utils.myPsOutput.colorTable[utils.myPsOutput.color[str(ci+1)]]
					for cc in c:
						print >> fa, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (cc, int(255*r),int(255*g),int(255*b))
				for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
					if (i1,i2) in tmpAretesApparent:
						print >> fa, "%s -- %s [style=\"dotted\"]" % (x[i1], x[i2])
					elif (i1,i2) in tmpAretesMemeEspece:
						print >> fa, "%s -- %s [style=\"dashed\"]" % (x[i1], x[i2])
					elif (i1,i2) in tmpAretes:
						print >> fa, "%s -- %s" % (x[i1], x[i2])

				print >> fa, "}"
				fa.close()


		nbA += len(clusters)


		# E. Ecriture
		for c in clusters:
			res.update(c)
			print >> f, " ".join(c)
			
			score = dict( [(e,[]) for e in geneBank.dicEspeces] )
			for i in c:
				score[geneBank.dicGenes[i][0]].append(i)
			l = [score[e][0] for e in score if len(score[e]) == 1]
			if len(l) >= 1:
				print >> ff, " ".join(l)
				nbO += 1
	
	# 3. On rajoute les genes qui n'ont plus qu'une seule copie
	for e in esp:
		for g in geneBank.dicEspeces[e].dicGenes:
			if (g not in lastComb) or (g in res):
				continue
			print >> f, g
			nbA += 1
			print >> ff, g
			nbO += 1
			res.add(g)

	f.close()
	ff.close()
	
	print >> sys.stderr, "OK (%d/%d)" % (nbA, nbO)
	
	return
	for (esp,_) in phylTree.items[anc]:
		if esp not in geneBank.dicEspeces:
			buildAncFile(esp, res)

buildAncFile(phylTree.root, set([]))
