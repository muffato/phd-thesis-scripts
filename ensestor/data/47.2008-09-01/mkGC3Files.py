#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import os
import sys
import time
import zipfile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs(  [("phylTree",file), ("GC3file",file)], [("ancGenesFiles",str,"")], "")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])

ancGenes = {}
for anc in phylTree.listAncestr:
	print >> sys.stderr, "Loading %s ..." % anc,
	genome = utils.myTools.myOpenFile(arguments["ancGenesFiles"] % phylTree.fileName[anc], "r")
	for gene in genome:
		(i,_,esp) = gene.replace('\n','').split('\t')
		ancGenes[frozenset(esp.split())] = (anc,i)
	genome.close()
	print >> sys.stderr, "OK"

f = utils.myTools.myOpenFile(arguments["GC3file"], "r")
for l in f:
	t = l.replace('\n','').split('\t')
	esp = frozenset(eval(t[1]))
	if esp in ancGenes:
		(anc,i) = ancGenes[esp]
		print "%s\t%s\t%d\t%s" % (anc,i,int(t[2])*3,t[3])
	else:
		print >> sys.stderr, esp
f.close()

