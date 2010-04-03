#! /users/ldog/muffato/python

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("CARs",file), ("humanGenome",file), ("identifier",str)], [], "")

humanGenome = utils.myGenomes.Genome(arguments["humanGenome"])

CARs = collections.defaultdict(list)
maHumanGenome = collections.defaultdict(list)

f = utils.myFile.openFile(arguments["CARs"], "r")
for l in f:
	l[:-1]
	if l.startswith("#"):
		CAR = l[1:-1]
	elif l.startswith(arguments["identifier"]):
		x = l[:-1].split()
		loc = x[0]
		strand = x[1]
		i = int(x[2][1:-1])
		x = loc.split(":")
		chrom = utils.myGenomes.commonChrName(x[0][8:].replace("chr",""))
		pos = tuple(int(x) for x in x[1].split("-"))
		assert pos[0] < pos[1]
		#for (j,gene) in enumerate(humanGenome.getGenesAt(chrom, pos[0], pos[1])):
		#	print utils.myFile.myTSV.printLine([i, j, j, gene.strand, gene.names[0]])
		#print i, strand, chrom, pos, " ".join([gene.names[0] for gene in humanGenome.getGenesAt(chrom, pos[0], pos[1], onlyInside=True)])
		maHumanGenome[chrom].append((pos,i))
		CARs[CAR].append((i,strand,chrom,pos))
f.close()

#sys.exit(0)

dicPos = {}
for chrom in maHumanGenome:
	maHumanGenome[chrom].sort()
	for (i,(_,j)) in enumerate(maHumanGenome[chrom]):
		dicPos[j] = i

dicStrand = {}
for g in humanGenome:
	dicStrand[g.names[0]] = g.strand

newCARs = {}
for CAR in CARs:
	#print CAR
	#print len(CARs[CAR])
	#print CARs[CAR]
	prev = CARs[CAR].pop(0)
	curr = list(prev[1:3] + prev[3])
	newCARs[CAR] = [curr]
	while len(CARs[CAR]) > 0:
		next = CARs[CAR].pop(0)
		#print prev, curr, next,
		if (next[2] == prev[2]) and (next[1] == prev[1]) and (((dicPos[next[0]] == dicPos[prev[0]]+1) and (prev[1] == "+")) or ((dicPos[next[0]] == dicPos[prev[0]]-1) and (prev[1] == "-"))):
			#print "extend", next[0], prev[0]
			#print "extend", prev[0], next[0]
			if prev[1] == "+":
				curr[3] = next[3][1]
			else:
				curr[2] = next[3][0]
		else:
			#print "stop"
			curr = list(next[1:3] + next[3])
			newCARs[CAR].append(curr)
		prev = next
	#print len(newCARs[CAR])
	#print newCARs[CAR]
	#print
	i = 0
	for x in newCARs[CAR]:
		#print CAR, utils.myFile.myTSV.printLine(x, delim=" "), " ".join([gene.names[0] for gene in humanGenome.getGenesAt(x[1], x[2], x[3], onlyInside=True)])
		#continue

		l = list(humanGenome.getGenesAt(x[1], x[2], x[3], onlyInside=True))
		if x[0] == "-":
			l.reverse()
		for g in l:
			print utils.myFile.myTSV.printLine([CAR, i, i, int(x[0]+"1")*dicStrand[g.names[0]], g.names[0]])
			i += 1
	



