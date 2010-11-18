#!/usr/bin/env python2

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("CARs",file), ("humanGenome",file), ("identifier",str)], [], "")

humanGenome = utils.myGenomes.Genome(arguments["humanGenome"])

#CARs = collections.defaultdict(list)
#maHumanGenome = collections.defaultdict(list)

f = utils.myFile.openFile(arguments["CARs"], "r")
for l in f:
	l[:-1]
	if l.startswith("#"):
		CAR = l[1:-1]
	elif l.startswith(arguments["identifier"]):
		x = l[:-1].split()
		loc = x[0]
		#strand = x[1]
		#i = int(x[2][1:-1])
		x = loc.split(":")
		chrom = utils.myGenomes.commonChrName(x[0][8:].replace("chr",""))
		pos = tuple(int(x) for x in x[1].split("-"))
		assert pos[0] < pos[1]
		#for (j,gene) in enumerate(humanGenome.getGenesAt(chrom, pos[0], pos[1])):
		#	print utils.myFile.myTSV.printLine([i, j, j, gene.strand, gene.names[0]])
		print chrom, pos, " ".join([gene.names[0] for gene in humanGenome.getGenesAt(chrom, pos[0], pos[1], onlyInside=True)])
		#maHumanGenome[chrom].append((pos,i))
		#CARs[CAR].append((i,strand,chrom,pos))
f.close()

