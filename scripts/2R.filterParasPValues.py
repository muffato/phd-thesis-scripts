#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit la matrice des p-values et extrait les clusters de fortement lies
"""

import sys
import bisect
import utils.myTools
import utils.myGenomes

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "paralogsList", "pvalueFile"], \
	[("cutoff",float,5)], \
	__doc__
)

f = utils.myTools.myOpenFile(noms_fichiers["pvalueFile"], "r")
pvalues = utils.myTools.defaultdict(dict)
for ligne in f:
	t = ligne.split()
	if t[0] == '+':
		chrom = t[1:]
	else:
		for (i,x) in enumerate(t[1:]):
			pvalues[t[0]][chrom[i]] = float(x)
f.close()

genome = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
paralogues = utils.myGenomes.Genome(noms_fichiers["paralogsList"])

print >> sys.stderr, "Filtering %d paralogs ..." % len(paralogues.lstGenes[None]),
nbPara = 0
for g in paralogues:
	tg = [c for (c,_) in genome.getPosition(g.names) if type(c) == int]
	chrpairs = frozenset(utils.myTools.myIterator.tupleOnStrictUpperList(tg))
	for (c1,c2) in chrpairs:
		if pvalues[str(c1)].get(str(c2),0) >= options["cutoff"]:
			print " ".join(g.names)
			break
print >> sys.stderr, "OK"

