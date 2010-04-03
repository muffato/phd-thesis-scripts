#! /users/ldog/muffato/python

__doc__ = """
	Separe les diagonales pairwise en 3 groupes en fonction de leur utilisation dans les blocs integres:
		1) celles utilisees a tous les ancetres intermediaires (= parfaites)
		2) celles non utilisees partout a cause de celles de 1) (= bloquees)
		3) le reste = celles non utilisees partout, mais a cause de 2) (= potentielles)
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import myDiags

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("target",str), ("pairwiseDiags",file),	("ancGenomes",str)], [], __doc__ )


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les ancetres cibles
targets = set(myDiags.getTargets(phylTree, arguments["target"])[1])

# Chargement des blocs integres
ancGenomes = {}
for anc in targets:
	ancGenomes[anc] = utils.myGenomes.Genome(arguments["ancGenomes"] % phylTree.fileName[anc])

# Est-ce que la paire est vue dans l'ancetre
def isAncOK(anc, (g1,s1), (g2,s2)):
	(c1,i1) = ancGenomes[anc].dicGenes[g1]
	(c2,i2) = ancGenomes[anc].dicGenes[g2]
	if c1 != c2:
		return False
	if i2 != (i1 + s1*ancGenomes[anc].lstGenes[c1][i1].strand):
		return False
	if s1*ancGenomes[anc].lstGenes[c1][i1].strand != s2*ancGenomes[anc].lstGenes[c2][i2].strand:
		return False
	return True


pairwiseAll = set()
pairwiseGood = set()
seen = set()

# Eclate un bloc pairwise en paires de genes et les teste pour chaque ancetre
def checkPerfectPairwise(t):
	diagS = [int(x) for x in t[7].split()]
	lanc = [x.split("=") for x in t[8].split("|")]
	lanc = []
	ldiagA = []
	for x in t[8].split("|"):
		(anc,diagA) = x.split("=")
		if anc in targets:
			lanc.append(anc)
			ldiagA.append(diagA.split())
	for ((lg1,s1,ea1,eb1),(lg2,s2,ea2,eb2)) in utils.myTools.myIterator.slidingTuple(zip(zip(*ldiagA), diagS, t[3].split(), t[6].split())):
		data = zip(lanc,lg1,lg2)
		t[0] = "2"
		t[3] = "%s %s" % (ea1,ea2)
		t[6] = "%s %s" % (eb1,eb2)
		t[7] = "%d %d" % (s1,s2)
		t[8] = "|".join( "%s=%s %s" % x for x in data )
		tt = tuple(t)
		pairwiseAll.add(tt)
		for (anc,g1,g2) in data:
			if not isAncOK(anc, (g1,s1), (g2,s2)):
				break
		else:
			print "\t".join(["PERFECT"] + t)
			pairwiseGood.add(tt)
			for (anc,g1,g2) in data:
				seen.add( (anc,g1,s1) )
				seen.add( (anc,g2,-s2) )

def checkBlockedPairwise(t):
	diagS = [int(x) for x in t[7].split()]
	lanc = [x.split("=") for x in t[8].split("|")]
	lanc = []
	ldiagA = []
	for x in t[8].split("|"):
		(anc,diagA) = x.split("=")
		if anc in targets:
			lanc.append(anc)
			ldiagA.append(diagA.split())
	((lg1,s1,ea1,eb1),(lg2,s2,ea2,eb2)) = zip(zip(*ldiagA), diagS, t[3].split(), t[6].split())
	for (anc,g1,g2) in zip(lanc,lg1,lg2):
		if (anc,g1,s1) in seen:
			break
		if (anc,g2,-s2) in seen:
			break
	else:
		return False
	return True


print >> sys.stderr, "Chargement des diagonales pairwise ...",
f = utils.myFile.myTSV.reader(arguments["pairwiseDiags"])
for t in f.csvobject:
	assert len(t) == 9
	checkPerfectPairwise(t)
f.file.close()
print >> sys.stderr, "%d/%d OK" % (len(pairwiseGood),len(pairwiseAll))

for t in pairwiseAll - pairwiseGood:
	if checkBlockedPairwise(t):
		print "\t".join(("BLOCKED",) + t)
	else:
		print "\t".join(("TOADD",) + t)

