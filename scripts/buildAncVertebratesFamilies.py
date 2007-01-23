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
	("orthoFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("paraFile",str,"~/work/data/paralogs/paras.%s.list.bz2")], \
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
	
	#for e in esp:
	#	doLoad(options["paraFile"] % e)
	#	utils.myTools.stderr.write('*')
		
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "...",
	#f = utils.myTools.stdout
	#ff = utils.myTools.null
	nbA = 0
	nbO = 0
	res = set([])
	for x in comb:
	
		# A. Composition du gene ancestral suivant les familles
		score = dict( [(e,[]) for e in geneBank.dicEspeces] )
		for g in x:
			if g not in geneBank.dicGenes:
				print >> sys.stderr, "GENE NON RECONNU: %s" % g,
				continue
			(e,_,_) = geneBank.dicGenes[g]
			score[e].append(g)
			
		# B. On filtre les genes qui ne sont pas specifiques a une branche
		if x[0] not in lastComb:
			nbBranchesOK = 0
			for (fils,_) in phylTree.items[anc]:
				if sum([len(score[e]) for e in phylTree.species[fils]]) >= 1:
					nbBranchesOK += 1
			if nbBranchesOK == 1:
				continue
		res.update(x)
	
		def test(i1, i2):
			
			g1 = x[i1]
			g2 = x[i2]

			if (g2 not in aretes[g1]) and (g1 not in aretes[g2]) and (geneBank.dicGenes[g1][0] != geneBank.dicGenes[g2][0]):
				return 0
			else:
				return 1

		#if len(x) > 1.5*len(esp)
		tmp = score.values()
		if (max([len(e) for e in tmp]) >= 3) or (tmp.count(2) > tmp.count(1)):
			#clusters = utils.myCommunities.launchCommunitiesBuild(len(x), test, minCoverage=0.9, minRelevance=0.4)
			clusters = utils.myCommunities.launchCommunitiesBuild(len(x), test, minCoverage=0.75, minRelevance=0.3)
		else:
			clusters = [range(len(x))]
		#tmp = utils.myCommunities.launchCommunitiesBuild(len(x), test, minCoverage=0.75, minRelevance=0.3)
		#clusters = utils.myCommunities.launchCommunitiesBuild(len(x), test)
		
		# ++ On verifie que le graphe est complet, sinon on l'affiche
		#complet = True
		#for (i,j) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		#	g1 = x[i]
		#	g2 = x[j]
		#
		#	if (g2 not in aretes[g1]) and (g1 not in aretes[g2]) and (geneBank.dicGenes[g1][0] != geneBank.dicGenes[g2][0]):
		#		complet = False
		#		break
			
		#for c in clusters:
		#	print " ".join([x[i] for i in c])
		#continue
		
		nbA += len(clusters)

		# C. Ecriture du gene ancestral
		#if len(clusters) > 1:
		if (sum([len(c) for c in clusters]) < len(x)) and (len(clusters) > 1):
		#if len(clusters) == 1 and len(tmp) > 1:
			f = open('/users/ldog/muffato/work/temp/graph/graph.%d' % nbA, 'w')
			print >> f, "graph {"
			for ci in xrange(len(clusters)):
				c = clusters[ci]
				(r,g,b) = utils.myPsOutput.colorTable[utils.myPsOutput.color[str(ci+1)]]
				for cc in c:
					print >> f, "%s [style=\"filled\",color=\"#%02X%02X%02X\"]" % (x[cc], int(255*r),int(255*g),int(255*b))
				
			for (i,j) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
				g1 = x[i]
				g2 = x[j]
	
				if (g2 in aretes[g1]) or (g1 in aretes[g2]) or (geneBank.dicGenes[g1][0] == geneBank.dicGenes[g2][0]):
					print >> f, "%s -- %s" % (g1, g2)
			print >> f, "}"
			f.close()
			#print >> f, " ".join(x)
			nbA += 1

		# D. Filtre des one2one
		#l = [score[e][0] for e in score if len(score[e]) == 1]
		#if len(l) >= 1:
		#	print >> ff, " ".join(l)
		#	nbO += 1
	
	# 3. On rajoute les genes qui n'ont plus qu'une seule copie
	#for e in esp:
	#	for g in geneBank.dicEspeces[e].dicGenes:
	#		if (g not in lastComb) or (g in res):
	#			continue
	#		print >> f, g
	#		nbA += 1
	#		print >> ff, g
	#		nbO += 1
	#		res.add(g)

	#f.close()
	#ff.close()
	
	print >> sys.stderr, "OK (%d/%d)" % (nbA, nbO)
	#del comb
	
	#for (esp,_) in phylTree.items[anc]:
	#	if esp not in geneBank.dicEspeces:
	#		buildAncFile(esp, res)

buildAncFile(phylTree.root, set([]))
