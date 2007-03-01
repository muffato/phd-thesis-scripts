#! /users/ldog/muffato/python -OO

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
import utils.myMaths
import utils.myTools


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"],\
	[("breakWhenFamilyNotComplete",bool,False), ("speciesList",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2")], \
	__doc__ \
)

esp = options["speciesList"].split(',')
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
#phylTree.loadSpeciesFromList(esp, options["genesFile"])
phylTree.loadAllSpeciesSince("Euteleostomi", options["genesFile"])


for l in sys.stdin:
	c = l.split()

	score = dict( [(e,0) for e in phylTree.dicGenomes] )

	for g in c:
		if g not in phylTree.dicGenes:
			if options["breakWhenFamilyNotComplete"]:
				break
			else:
				continue
		(e,_,_) = phylTree.dicGenes[g]
		score[e] += 1
	else:
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

