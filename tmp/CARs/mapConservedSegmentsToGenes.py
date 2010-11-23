#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("ConservedSegments",file), ("CARs",file), ("humanGenome",file), ("dogGenome",file), ("boreoGenome",file)], \
	[("onlyInside",bool,True)], \
	"" \
)

humanGenome = utils.myGenomes.Genome(arguments["humanGenome"])
dogGenome = utils.myGenomes.Genome(arguments["dogGenome"])
boreoGenome = utils.myGenomes.Genome(arguments["boreoGenome"])

def getGenesFrom(genome, pos):
	l = list(genome.getGenesAt(pos[0], pos[1][0], pos[1][1], onlyInside=arguments["onlyInside"]))
	if pos[2] == "+":
		return [(g.names[0],g.strand) for g in l]
	else:
		return [(g.names[0],-g.strand) for g in reversed(l)]

def findFirstCommonGene(lH, lD):
	all = []
	for (iH,(_,posH)) in enumerate(lH):
		for (iD,(_,posD)) in enumerate(lD):
			if posH == posD:
				all.append((iH+iD, iH, iD))
	return min(all) if len(all) > 0 else None

def do(segID, dataH, dataD):
	lH = getGenesFrom(humanGenome, dataH)
	lHb = [(g,boreoGenome.dicGenes[g[0]]) for g in lH if g[0] in boreoGenome.dicGenes]
	lD = getGenesFrom(dogGenome, dataD)
	lDb = [(g,boreoGenome.dicGenes[g[0]]) for g in lD if g[0] in boreoGenome.dicGenes]
	
	res = [segID]
	
	#print segID
	#print "human", dataH
	#print len(lH), lH
	#print len(lHb), lHb
	#print "dog", dataD
	#print len(lD), lD
	#print len(lDb), lDb
	
	i1 = findFirstCommonGene(lHb, lDb)
	if i1 is not None:
		res.extend(lHb[i1[1]][0])
		res.extend(lDb[i1[2]][0])
		r1 = (lHb[i1[1]][0], lDb[i1[2]][0])
		#print "common1", i1, lH[i1[1]], lD[i1[2]]
	else:
		res.extend([None]*4)
		r1 = None
		#print "common1", None

	lHb.reverse()
	lDb.reverse()
	i2 = findFirstCommonGene(lHb, lDb)
	#i2 = findFirstCommonGene(list(reversed(lH)), list(reversed(lD)))
	if i2 is not None:
		res.extend(lHb[i2[1]][0])
		res.extend(lDb[i2[2]][0])
		r2 = (lHb[i2[1]][0], lDb[i2[2]][0])
		#print "common2", i2, lH[i2[1]], lD[i2[2]]
	else:
		res.extend([None]*4)
		r2 = None
		#print "common2", None

	#print utils.myFile.myTSV.printLine(res)
	return (r1, r2)


def readCoords(s):
	t = s.split()
	(_,_,s) = t[0].partition(".")
	(chrom,s) = s.split(":")
	chrom = utils.myGenomes.commonChrName(chrom[3:])
	return (chrom, tuple(int(x) for x in s.split("-")), t[1])

mapping = {}

f = utils.myFile.openFile(arguments["ConservedSegments"], "r")
segID = None
for l in f:
	if l.startswith(">"):
		if segID is not None:
			mapping[segID] = do(segID, dataH, dataD)
		segID = int(l[1:-1])
	elif l.startswith("hg18"):
		dataH = readCoords(l)
	elif l.startswith("canFam2"):
		dataD = readCoords(l)
if segID is not None:
	mapping[segID] = do(segID, dataH, dataD)
f.close()

#sys.exit(0)

def revMapping(mi):
	return ((mi[0][0],-mi[0][1]), (mi[1][0],-mi[1][1]))

seen = set()
CARs = collections.defaultdict(list)
f = utils.myFile.openFile(arguments["CARs"], "r")
for l in f:
	l[:-1]
	if l.startswith("#"):
		CAR = l[1:-1]
	elif l.startswith("hg18"):
		x = l[:-1].split()
		i = int(x[2][1:-1])
		strand = x[1]
		
		if (mapping[i][0] is not None) and (mapping[i][1] is not None):
			if strand == "-":
				mapping[i] = (revMapping(mapping[i][1]), revMapping(mapping[i][0]))
			CARs[CAR].append(i)
			seen.add(mapping[i][0][0][0])
			seen.add(mapping[i][0][1][0])
			seen.add(mapping[i][1][0][0])
			seen.add(mapping[i][1][1][0])
f.close()
print >> sys.stderr, sum(len(x) for x in CARs.itervalues()), len(seen)

def rewriteGenome(genome, seen):
	dicPos = {}
	dicLen = {}
	for chrom in genome.lstGenes:
		i = 0
		for g in genome.lstGenes[chrom]:
			# Uniquement les genes presents dans un segment conserve
			if len(seen.intersection(g.names)) > 0:
				for s in seen.intersection(g.names):
					dicPos[s] = (chrom, i, g.strand)
				i += 1
		dicLen[chrom] = i
	return (dicPos, dicLen)


(humanPos, humanLen) = rewriteGenome(humanGenome, seen)
(dogPos, dogLen) = rewriteGenome(dogGenome, seen)
(boreoPos, boreoLen) = rewriteGenome(boreoGenome, seen)

print >> sys.stderr, len(humanPos), len(dogPos), len(boreoPos)

def getStatus((g1,s1), (g2,s2), dicPos, dicLen, genome):
	(c1,i1,t1) = dicPos[g1]
	(c2,i2,t2) = dicPos[g2]
	l1 = genome.getOtherNames(g1)
	l2 = genome.getOtherNames(g2)
	print dicPos[g1], dicPos[g2], g1 if len(l1) == 0 else l1[0], g2 if len(l2) == 0 else l2[0], dicLen[c1], dicLen[c2]

	if c1 != c2:
		if ((s1 == t1) and (i1 == dicLen[c1]-1)) or (i1 == 0):
			if ((s2 == t2) and (i2 == 0)) or (i2 == dicLen[c2]-1):
				# 2 extremites de nos blocs
				return 1
	else:
		if (s1 == t1) and (s2 == t2) and (i2 == i1+1):
			# Adjacence presente chez nous
			return 2
		if (s1 == -t1) and (s2 == -t2) and (i2 == i1-1):
			# Adjacence presente chez nous
			return 2
		if i1 == i2:
			return -1

	# Difference
	return 0

for CAR in CARs:
	for (seg1,seg2) in utils.myTools.myIterator.slidingTuple(CARs[CAR]):

		print
		print seg1, seg2
		sH = getStatus(mapping[seg1][1][0], mapping[seg2][0][0], humanPos, humanLen, humanGenome)
		sD = getStatus(mapping[seg1][1][1], mapping[seg2][0][1], dogPos, dogLen, dogGenome)
		sA = getStatus(mapping[seg1][1][0], mapping[seg2][0][0], boreoPos, boreoLen, boreoGenome)
		print "hum", sH
		print "dog", sD
		print "anc", sA
		print "tot", "hum%d" % sH, "dog%d" % sD, "anc%d" % sA, seg1, seg2
		if ((sD == 1) or (sD == 0)) and (sA == 2):
			print "rearr-dog", boreoGenome.getOtherNames(mapping[seg1][1][1][0])[0], boreoGenome.getOtherNames(mapping[seg2][0][1][0])[0]
		if ((sH == 1) or (sH == 0)) and (sA == 2):
			print "rearr-hum", boreoGenome.getOtherNames(mapping[seg1][1][0][0])[0], boreoGenome.getOtherNames(mapping[seg2][0][0][0])[0]

