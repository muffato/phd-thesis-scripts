#! /users/ldog/muffato/python

__doc__ = """
Lit des familles de genes sur l'entree standard.
Calcule le nombre de genes de chaque espece et agit en consequence
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPhylTree


########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)],\
	[("breakWhenFamilyNotComplete",bool,False), ("speciesList",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
if arguments["speciesList"] == "":
	phylTree.loadAllSpeciesSince(None, arguments["genesFile"])
else:
	phylTree.loadSpeciesFromList(arguments["speciesList"].split(','), arguments["genesFile"])


for e in phylTree.dicGenomes:
	print "%s\t" % e,
print

for l in sys.stdin:
	c = l.split()

	score = dict( [(e,0) for e in phylTree.dicGenomes] )

	for g in c:
		if g not in phylTree.dicGenes:
			if arguments["breakWhenFamilyNotComplete"]:
				break
			else:
				print >> sys.stderr, "Can't find %s" % g
				continue
		(e,_,_) = phylTree.dicGenes[g]
		score[e] += 1
	else:

		if min(score.values()) >= 1:
			print l,
		continue

		for e in phylTree.dicGenomes:
			print "%d\t" % score[e],
		print
		continue

		t = score.values()
		t = [x for x in score.values() if x != 0]
		t.sort()
		if t[len(t)/2] >= 10:
			print score
			#print l,
		#print sum(t)/(len(t)-t.count(0))
		#t = [score[x] for x in phylTree.lstEspecesNonDup]
		#tt = [score[x] for x in phylTree.lstEspecesDup]
		#tt = [score[x] for x in ]
		#if score['H'] == 1 and score['C'] == 1:
		#if max(t) > 0 and max(tt) > 0:
		#if max(t) == 2 and max(tt) == 1:
		#if sum(t) < 13 and sum(tt) < 6:
		#if len(t) == 1 and 1 in t:
		#if max(t) <= 1: # and max(tt) <= 2:
		#if max(t) == 1 and min(t) == 1:
		#	print l,

