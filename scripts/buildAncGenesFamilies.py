#! /usr/bin/python2.4

__doc__ = """
Lit l'arbre phylogenetique et construit les fichiers de genes ancestraux
 pour tous les noeuds en utilisant les fichiers d'orthologues de Biomart
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("orthoFile",str,"/users/ldog/muffato/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("one2oneFile",str,"/users/ldog/muffato/work/data/ancGenes/one2one.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.getSpecies(phylTree.root))

for anc in phylTree.items:
	esp = phylTree.getSpecies(anc)
	ligne = "bzcat"
	for e1 in esp:
		for e2 in esp:
			if e1 == e2:
				continue
			ligne += " " + (options["orthoFile"] % (e1,e2))
	ligne += " | awk '{print $1, $4}' | /users/ldog/muffato/work/scripts/groupObjects.py | bzip2 > " + (options["ancGenesFile"] % anc)
	
	print "Familles d'orthologues de", anc, ":", esp
	os.system(ligne)
	
	print "Familles de one2one de", anc, ":", esp

	f = utils.myTools.myOpenFile(options["ancGenesFile"] % anc, 'r')
	ff = utils.myTools.myOpenFile(options["one2oneFile"] % anc, 'w')
	
	for l in f:
		c = l.split()

		score = dict( [(e,[]) for e in geneBank.dicEspeces] )
		
		for g in c:
			if g not in geneBank.dicGenes:
				print >> sys.stderr, "ERREUR:", g
				continue
			(e,_,_) = geneBank.dicGenes[g]
			score[e].append(g)
			
		l = [score[e][0] for e in score if len(score[e]) == 1]
		if len(l) >= 2:
			print >> ff, " ".join(l)
			
	f.close()
	ff.close()

