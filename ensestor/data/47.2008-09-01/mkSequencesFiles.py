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

arguments = utils.myTools.checkArgs(  [("phylTree",file)], [("range",str,""), ("ancGenesFiles",str,""),("ancSequenceFiles",str,"")], "")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])

ancGenes = collections.defaultdict(list)
for anc in phylTree.listAncestr:
	print >> sys.stderr, "Loading %s ..." % anc,
	genome = utils.myTools.myOpenFile(arguments["ancGenesFiles"] % phylTree.fileName[anc], "r")
	for gene in genome:
		(i,name,esp) = gene.replace('\n','').split('\t')
		ancGenes[frozenset(esp.split())].append( (anc,i,name) )
	genome.close()
	print >> sys.stderr, "OK"

for i in utils.myTools.getRange(arguments["range"]):
	print >> sys.stderr, i,
	f = utils.myTools.myOpenFile(arguments["ancSequenceFiles"] % i, "r")
	#f =  zipfile.ZipFile(arguments["IN.familyZipFile"] % i, "r")
	#for l in f.read("ancSeq3.%d.txt" % i).split('\n'):
	#for l in f.read("ancSequence.txt").split('\n'):
	for l in f:
		t = l.replace('\n','').split('\t')
		if t[0] == "SPECIES":
			anc = ancGenes[frozenset(t[1].split())]
			if anc == []:
				sys.stderr.write("!")
		elif (t[0] == "SEQ"): # and (anc != None):
			#seq = "NN" + "NN".join(t[1].replace('-',''))
			seq = t[1].replace('-','')
			for (anc,i,name) in anc:
				print "%s\t%s\t%s\t%s" % (anc,i,name,seq)
	f.close()
	print >> sys.stderr, " OK"

