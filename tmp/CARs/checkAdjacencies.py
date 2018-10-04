#!/usr/bin/env python2

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("CARs",file), ("humanGenome",file), ("boreoGenome",file), ("identifier",str)], [], "")

humanGenome = utils.myGenomes.Genome(arguments["humanGenome"])
boreoGenome = utils.myGenomes.Genome(arguments["boreoGenome"])

seenHuman = set()

CARs = collections.defaultdict(list)

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

		# La liste des genes humains dans l'intervalle ...
		l = list(humanGenome.getGenesAt(chrom, pos[0], pos[1]))
		# ... totalement inclus ...
		#l = [g for g in l if g.beginning >= pos[0] and g.end <= pos[1]]
		# ... et presents chez l'ancetre
		l = [g for g in l if g.names[0] in boreoGenome.dicGenes]
		# ... et pas en singletons
		#l = [g for g in l if len(boreoGenome.lstGenes[boreoGenome.dicGenes[g.names[0]][0]]) >= 2]


		if len(l) > 0:

			# Pour reduire les genomes humains/boreo
			seenHuman.update(g.names[0] for g in l)

			# Les bornes dans le bon sens
			if strand == "+":
				l = [(g.names[0],g.strand) for g in l]
			else:
				l = [(g.names[0],-g.strand) for g in reversed(l)]
		
			CARs[CAR].append( (i, l[0], l[-1]) )
f.close()

print >> sys.stderr, sum(len(x) for x in CARs.itervalues())

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


(humanPos, _) = rewriteGenome(humanGenome, seenHuman)
(boreoPos, boreoLen) = rewriteGenome(boreoGenome, seenHuman)

print >> sys.stderr, len(seenHuman), len(humanPos), len(boreoPos)

def getStatus(g1, g2, dicPos):
	(c1,i1,t1) = dicPos[g1]
	(c2,i2,t2) = dicPos[g2]

	if c1 != c2:
		if ((s1 == t1) and (i1 == boreoLen[c1]-1)) or (i1 == 0):
			if ((s2 == t2) and (i2 == 0)) or (i2 == boreoLen[c2]-1):
				# 2 extremites de nos blocs
				return 1
	else:
		if (s1 == t1) and (s2 == t2) and (i2 == i1+1):
			# Adjacence presente chez nous
			return 2
		if (s1 == -t1) and (s2 == -t2) and (i2 == i1-1):
			# Adjacence presente chez nous
			return 2

	# Difference
	return 0


for CAR in CARs:
	for ((seg1,_,(g1,s1)),(seg2,(g2,s2),_)) in utils.myTools.myIterator.slidingTuple(CARs[CAR]):
		#assert g1 != g2, (seg1,seg2,(g1,s1),(g2,s2))
		if g1 == g2:
			continue
		(c1,i1,t1) = boreoPos[g1]
		(c2,i2,t2) = boreoPos[g2]

		print
		print (g1,s1), (g2,s2)
		print humanPos[g1], humanPos[g2]
		print boreoPos[g1], boreoPos[g2]
		print boreoLen[c1], boreoLen[c2]
		print humanGenome.dicGenes[g1], humanGenome.dicGenes[g2]
		print boreoGenome.dicGenes[g1], boreoGenome.dicGenes[g2]
		print "anc" + str(getStatus(g1, g2, boreoPos)), "hum" + str(getStatus(g1, g2, humanPos)), seg1, seg2

