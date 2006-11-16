#! /users/ldog/muffato/python

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


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("homologyLevels",str,"ortholog_one2many,ortholog_many2many,apparent_ortholog_one2one,ortholog_one2one"), \
	("orthoFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("one2oneFile",str,"~/work/data/one2one/one2one.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.getSpecies(phylTree.root))
homologies = options["homologyLevels"].split(",")

def buildAncFile(anc, lastComb):

	# 1. on combine tous les fichiers d'orthologues
	comb = utils.myTools.myCombinator([])
	esp = phylTree.getSpecies(anc)
	print >> sys.stderr, "Construction des familles d'orthologues de", anc, ":", "-".join(esp), "",
	for (i,j) in utils.myTools.myMatrixIterator(len(esp), len(esp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		f = utils.myTools.myOpenFile(options["orthoFile"] % (esp[i],esp[j]), 'r')
		for ligne in f:
			c = ligne.split()
			if c[6] not in homologies:
				continue
			comb.addLink([c[0], c[3]])
		f.close()
		sys.stderr.write('.')
	print >> sys.stderr, " OK"
	
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "...",
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
				print >> sys.stderr, "GENE NON RECONNU: %s" % g,
				continue
			(e,_,_) = geneBank.dicGenes[g]
			score[e].append(g)
			
		# B. On filtre les genes qui ne sont pas specifiques a une branche
		if x[0] not in lastComb:
			nbBranchesOK = 0
			for (fils,_) in phylTree.items[anc]:
				if sum([len(score[e]) for e in phylTree.getSpecies(fils)]) >= 1:
					nbBranchesOK += 1
			if nbBranchesOK == 1:
				continue
		res.update(x)
		
		# C. Ecriture du gene ancestral
		print >> f, " ".join(x)
		nbA += 1

		# D. Filtre des one2one
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
	del comb
	
	for (esp,_) in phylTree.items[anc]:
		if esp not in geneBank.dicEspeces:
			buildAncFile(esp, res)

buildAncFile(phylTree.root, set([]))
