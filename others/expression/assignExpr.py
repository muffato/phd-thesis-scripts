#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import collections

import utils.myMaths
import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("exprFile",file), ("blatFile",file)], [], "Cree une liste ordonnee des genes en tenant compte du plus petit transcrit")

dicGP = {}
f = utils.myFile.openFile(arguments["blatFile"], "r")
for l in f:
	t = l.split()
	lg = [x for x in t if x.startswith("ENSG0")]
	assert len(lg) == 1
	g = lg[0]
	assert g not in dicGP
	t2 = [x for x in t if not x.startswith("ENS")]
	if len(t2) != 0:
		dicGP[g] = t2
f.close()

dicExpr = {}
names = None
f = utils.myFile.openFile(arguments["exprFile"], "r")
for l in f:
	t = l.split()
	if names == None:
		assert len(t) == 2*int(len(t)/2)
		names = [t[2*i] for i in xrange(len(t)/2)]
	else:
		#dicExpr[t[0]] = [int(x) for x in t[1:]]
		dicExpr[t[0]] = t[1:]
f.close()

print "GENE_NAME", "PROBESET_NAME", " ".join(names)
for (g,lp) in dicGP.iteritems():
	lp = [p for p in lp if p in dicExpr]
	assert len(lp) > 0
	for p in lp:
		print g, p, " ".join([dicExpr[p][2*i] for i in xrange(len(names))])
		print g, p, " ".join([dicExpr[p][2*i+1] for i in xrange(len(names))])
	continue
	if len(lp) == 0:
		continue
	for i in xrange(len(names)):
		le = [dicExpr[p][2*i] for p in lp] + [dicExpr[p][2*i+1] for p in lp]
		s = utils.myMaths.myStats.valSummary(le)
		print g, i, s[2], s[8], s[9], s[0], s[7], le

