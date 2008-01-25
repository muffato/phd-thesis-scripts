#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree", "tree.5"], \
	[("IN.CDS",str,""), ("OUT.FASTA-CDS",str,""), ("OUT.FASTA-Prot",str,""), ("OUT.tree",str,""), ("transeqBin",str,"")], \
	"Ecrit les familles" \
)

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree"], buildLinks = False)

CDS = {}
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Chargement des CDS de", esp, "...",
	f = utils.myTools.myOpenFile(options["IN.CDS"] % phylTree.fileName[esp], "r")
	for ligne in f:
		t = ligne.split('\t')
		CDS[t[0]] = t[1]
	f.close()
	print >> sys.stderr, "OK"


f = utils.myTools.myOpenFile(noms_fichiers["tree.5"], "r")
n = 0
for (i,ligne) in enumerate(f):
	print >> sys.stderr, "Ligne", i,
	if (i % 2) == 0:
		n += 1
		utils.myTools.mkDir(options["OUT.FASTA-CDS"] % n)
		fasta = utils.myTools.myOpenFile(options["OUT.FASTA-CDS"] % n, "w")
		for g in ligne.split():
			print >> fasta, ">" + g
			print >> fasta, CDS[g]
		fasta.close()
		os.system("%s -sequence %s -trim | nawk '$1~/>/{gsub(/_1/,\"\")}{print}' > %s" % (options["transeqBin"],options["OUT.FASTA-CDS"] % n,options["OUT.FASTA-Prot"] % n))
	else:
		tree = utils.myTools.myOpenFile(options["OUT.tree"] % n, "w")
		print >> tree, ligne,
		tree.close()
	print >> sys.stderr, "OK"
f.close()
