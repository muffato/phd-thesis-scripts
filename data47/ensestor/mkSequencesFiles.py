#! /users/ldog/muffato/python -OO

import os
import sys
import time
import zipfile
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs(  ["phylTree"], [("start",int,0), ("end",int,0), ("ancGenesFiles",str,""),("IN.familyZipFile",str,"")], "")

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

for i in xrange(options["start"], options["end"]+1):
	print >> sys.stderr, i,
	f =  zipfile.ZipFile(options["IN.familyZipFile"] % i, "r")
	#for l in f.read("ancSeq3.%d.txt" % i).split('\n'):
	for l in f.read("ancSequence.txt").split('\n'):
		t = l.split('\t')
		if t[0] == "SPECIES":
			esp = frozenset(t[1].split())
			if esp in ancGenes:
				(anc,i) = ancGenes[esp]
			else:
				print >> sys.stderr, t[1], "NOT FOUND",
				anc = None
		elif (t[0] == "SEQ") and (anc != None):
			seq = "NN" + "NN".join(t[1].replace('-',''))
			print "%s\t%s\t%s" % (anc,i,seq)
	f.close()
	print >> sys.stderr, "OK"

