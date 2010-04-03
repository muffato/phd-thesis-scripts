#! /users/ldog/muffato/python

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("range",str), ("dicGeneProtTrans",str)], \
	[("ancGenesFile",str,""), ("mattFile",str,""), ("outputFile",str,"")], \
	"Extrait les sequences de chaque ancetre" \
)


# Association gene_name > transcript_name
f = utils.myFile.openFile(arguments["dicGeneProtTrans"], "r")
dicProtTransc = {}
for l in f:
	t = l.split()
	dicProtTransc[t[0]] = t[2]
f.close()

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

lancGenes = {}
files = {}
for e in phylTree.listAncestr:
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[e])
	for g in ancGenes:
		s = frozenset(dicProtTransc[x] for x in g.names[1:])
		if s in lancGenes:
			if e in phylTree.allDescendants[lancGenes[s][0]]:
				lancGenes[s] = (e, g.names[0])
		else:
			lancGenes[s] = (e, g.names[0])
	files[e] = utils.myFile.openFile(arguments["outputFile"] % phylTree.fileName[e], "w")

for i in utils.myTools.getRange(arguments["range"]):
	try:
		rec = utils.myFile.openFile(arguments["mattFile"] % i, "r")
	except IOError:
		print >> sys.stderr, "NOT FOUND:", arguments["mattFile"] % i
		continue
	for data in rec:
		data = data[:-1]
		if data.startswith("SPECIES"):
			inispec = data.split("\t")[1]
			names = frozenset(inispec.split())
		elif data.startswith("SEQ"):
			#print i,
			if len(names) == 1:
				#print "SINGLE", names
				continue
			assert len(names) >= 2
			if names not in lancGenes:
				#print "NOTFOUND", names
				continue
			#print "FOUND", names
			print utils.myFile.myTSV.printLine([i, lancGenes[names][0], lancGenes[names][1], inispec])
			utils.myGenomes.printFastaSeq(files[lancGenes[names][0]], "-".join(lancGenes[names]), data.split("\t")[1].replace('-', ''))

for e in phylTree.listAncestr:
	files[e].close()

