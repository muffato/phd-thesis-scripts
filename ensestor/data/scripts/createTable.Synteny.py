#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Cree la base Genomicus (table Synteny)
	Necessite la table Gene pour positionner les blocs sur les chromosomes
	Lit les fichiers de diagonales + genes ancestraux
"""

import sys
import collections

import utils.myTools
import utils.myGenomes
import utils.myPhylTree



def storeAncDiags(anc):

	ancGenes = utils.myGenomes.Genome(arguments["ancGenes"] % phylTree.fileName[anc])
	genome = utils.myGenomes.Genome(arguments["diags"] % phylTree.fileName[anc], ancGenes=ancGenes)
	global synt_id

	# Parcours des chromosomes
	print >> sys.stderr, "Inserting ancestral diags", anc, "...",
	for (chrom,l) in genome.lstGenes.iteritems():
		if len(l) == 1:
			continue
		trans = [dicPos[phylTree.indNames[anc]][gene.names[0]] for gene in l]
		transS = sorted(trans)
		# Tout le monde est sur un seul chromosome
		assert transS[0][0] == transS[-1][0]
		assert transS[0][1] + len(trans) - 1 == transS[-1][1]
		if trans != transS:
			assert list(transS.__reversed) == transS
		print utils.myFile.myTSV.MySQLFileWriter((synt_id, phylTree.indNames[anc], trans[0][0], transS[0][1], transS[-1][1]))
		synt_id += 1
	print >> sys.stderr, "OK"


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("GeneTable",file)], \
	[("ancGenes",str,""), ("diags",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Lecture de la table des genes
print >> sys.stderr, "Loading table 'Gene' ...",
dicPos = collections.defaultdict(dict)
for line in utils.myFile.myTSV.readTabular(arguments["GeneTable"], [str]*10):
	dicPos[int(line[1])][line[2]] = (line[3],int(line[4]))
print >> sys.stderr, "OK"

synt_id = 0
for anc in phylTree.listAncestr:
	storeAncDiags(anc)



