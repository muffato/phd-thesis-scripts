#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des arbres de proteines
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("ensemblTree",file)], [("cdsFile",str,"")], __doc__)

dicTr = {}
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
for e in phylTree.listSpecies:
	f = arguments["cdsFile"] % e.replace(' ', '_')
	if utils.myFile.hasAccess(f):
		print >> sys.stderr, e
		dicTr.update((n.split()[0],s) for (n,s) in utils.myGenomes.myFASTA.iterFile(f))
		#dicTr.update(utils.myGenomes.myFASTA.iterFile(f))

print >> sys.stderr, "OK"
for (i,(root,data,info)) in enumerate(utils.myProteinTree.loadTree(arguments["ensemblTree"])):
	for d in info.itervalues():
		if 'protein_alignment' in d:
			#print i, d['gene_name'], d['protein_name'], d['transcript_name'], d['protein_alignment'], dicTr["%(gene_name)s|%(transcript_name)s" % d]
			print i, d['gene_name'], d['protein_name'], d['transcript_name'], d['protein_alignment'], dicTr[d['protein_name']]


