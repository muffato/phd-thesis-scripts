#! /users/ldog/muffato/python -OO

import os
import sys
import time
import zipfile
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs(  ["phylTree", "GC3file"], [("ancGenesFiles",str,"")], "")

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree"])

ancGenes = {}
for anc in phylTree.listAncestr:
	print >> sys.stderr, "Loading %s ..." % anc,
	genome = utils.myTools.myOpenFile(options["ancGenesFiles"] % phylTree.fileName[anc], "r")
	for gene in genome:
		(i,_,esp) = gene[:-1].split('\t')
		ancGenes[frozenset(esp.split())] = (anc,i)
	genome.close()
	print >> sys.stderr, "OK"

f = utils.myTools.myOpenFile(noms_fichiers["GC3file"], "r")
for l in f:
	t = l[:-1].split('\t')
	esp = frozenset(eval(t[1]))
	if esp in ancGenes:
		(anc,i) = ancGenes[esp]
		print "%s\t%s\t%d\t%s" % (anc,i,int(t[2])*3,t[3])
	else:
		print >> sys.stderr, esp
f.close()

