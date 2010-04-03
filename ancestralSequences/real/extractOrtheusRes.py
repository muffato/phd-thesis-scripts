#! /users/ldog/muffato/python

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("range",str)], \
	[("genesFile",str,""), ("ancGenesFile",str,""), ("ortheusFile",str,""), ("outputFile",str,"")], \
	"Extrait les sequences de chaque ancetre" \
)


# Association gene_name > transcript_name
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
notransc = {}
for e in phylTree.listSpecies:
	f = utils.myFile.openFile(arguments["genesFile"] % phylTree.fileName[e], "r")
	for l in f:
		t = l.split("\t")
		if len(t) >= 7:
			notransc[t[4]] = t[6]
	f.close()

lancGenes = {}
for e in phylTree.listAncestr:
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[e])
	for g in ancGenes:
		s = frozenset(notransc[x] if x in notransc else x for x in g.names[1:])
		if s in lancGenes:
			if e in phylTree.allDescendants[lancGenes[s][0]]:
				lancGenes[s] = (e, g.names[0])
		else:
			lancGenes[s] = (e, g.names[0])

for i in utils.myTools.getRange(arguments["range"]):
	try:
		fasta = utils.myGenomes.myFASTA.loadFile(arguments["ortheusFile"] % i)
	except IOError:
		print >> sys.stderr, "NOT FOUND:", arguments["ortheusFile"] % i
		continue
	for s in fasta:
			print i,
			names = frozenset(s.split('_'))
			if len(names) == 1:
				print "SINGLE", names
				continue
			assert len(names) >= 2
			if names not in lancGenes:
				print "NOTFOUND", names
				continue
			print "FOUND", names
			f = utils.myFile.openFile(arguments["outputFile"] % ((i,) + lancGenes[names]), "w")
			utils.myGenomes.printFastaSeq(f, "%s/%s" % lancGenes[names], fasta[s].replace('-', ''))
			f.close()


