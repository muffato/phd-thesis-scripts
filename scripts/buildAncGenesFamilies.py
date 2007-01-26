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
			comb.addLink([gA, gB])
			if gA not in aretes:
				aretes[gA] = [gB]
			else:
				aretes[gA].append(gB)
			if gB not in aretes:
				aretes[gB] = [gA]
			else:
				aretes[gB].append(gA)
		f.close()
		


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
		score = dict( [(e,[]) for e in geneBank.dicEspeces] )
		for g in x:
			if g not in geneBank.dicGenes:
				print >> sys.stderr, "GENE NON RECONNU: %s" % g
				continue
			(e,_,_) = geneBank.dicGenes[g]
			score[e].append(g)
			
		poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
		poidsTout = score.values()
		
		# B. On enleve les nouvelle familles qui ne sont pas dans toutes les sous-branches
		if (x[0] not in lastComb) and (0 in poidsBranches):
			continue

		# La fonction qui renvoie le s
		def test(i1, i2):
			if (i1,i2) in tmpAretes:
				return 1
			else:
				return 0

		# C. On calcule les aretes du graphe
		tmpAretes = []
		tmpAretesMemeEspece = []
		for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			g1 = x[i1]
			g2 = x[i2]
			if (g1 in aretes[g2]):
				tmpAretes.append( (i1,i2) )
			if (geneBank.dicGenes[g1][0] == geneBank.dicGenes[g2][0]):
				tmpAretesMemeEspece.append( (i1,i2) )


		# D. Filtres
		# D.1/ Graphe complet -> 1 cluster
		if (len(tmpAretes) + len(tmpAretesMemeEspece)) == (len(x)*(len(x)-1))/2:
			print >> sys.stderr, "Graphe complet"
			clusters = [x]
			
		# D.2/ Tous les genes en 1 copie dans une des branches > 1 cluster
		elif min(poidsBranches) == 1:
			print >> sys.stderr, "Branche sans duplication"
			clusters = [x]

		# D.3/ On est oblige de clusteriser
		else:
			#(relev,clusters) = utils.myCommunities.launchCommunitiesBuild(len(x), test, minCoverage=0.9, minRelevance=0.4)
			(relev,clusters) = utils.myCommunities.launchCommunitiesBuild(len(x), test)
			clusters = [[x[i] for i in c] for c in clusters]

			if relev[0] < 0.4 or len(utils.myMaths.flatten(clusters)) <= 0.9*len(x):
				fa = open('/users/ldog/muffato/work/tutu/graph-%f' % relev[0], 'w')
				print >> fa, "graph {"
				for ci in xrange(len(clusters)):
					c = clusters[ci]
					(r,g,b) = utils.myPsOutput.colorTable[utils.myPsOutput.color[str(ci+1)]]
					for cc in c:
						print >> fa, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (cc, int(255*r),int(255*g),int(255*b))
				for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
					if (i1,i2) in tmpAretes:
						print >> fa, "%s -- %s" % (x[i1], x[i2])
				print >> fa, "}"
				fa.close()


		nbA += len(clusters)

		# @. Ecriture du gene ancestral
		#if len(clusters) >= 2:
		#	f = open('/users/ldog/muffato/work/temp/graph/graph-%f' % relev[0], 'w')
		#	print >> f, "graph {"
		#	for ci in xrange(len(clusters)):
		#		c = clusters[ci]
		#		(r,g,b) = utils.myPsOutput.colorTable[utils.myPsOutput.color[str(ci+1)]]
		#		for cc in c:
		#			print >> f, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (x[cc], int(255*r),int(255*g),int(255*b))
		#	for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		#		if (i1,i2) in tmpAretes:
		#			print >> f, "%s -- %s" % (x[i1], x[i2])
		#	print >> f, "}"
		#	f.close()
		#	#print >> f, " ".join(x)

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
