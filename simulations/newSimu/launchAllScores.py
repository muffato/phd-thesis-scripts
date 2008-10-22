#! /users/ldog/muffato/python

import os
import sys
import itertools

import utils.myPhylTree

Executable = '/workspace/muffato/simulations/newSimu/computeScores.py'

phylTree = utils.myPhylTree.PhylogeneticTree("/workspace/muffato/data49/phylTree.noLowCoverage.49.conf")

anc = sys.argv[1].replace('^', ' ')
sign = sys.argv[2]

for (scoring,weightM,weightP) in itertools.product(range(5),"-+","-+"):
#for (scoring,weightM,weightP) in itertools.product(range(1),"-","-"):

	ancF = phylTree.fileName[anc]
	BaseRep = "/workspace2/muffato/newSimu/all2/signatures/"
	BaseName = BaseRep + "scores/%s.%s.%d.%s.%s" % (sign,ancF,scoring,weightM,weightP)
	Arguments = "+psyco +bz2 /workspace/muffato/data49/phylTree.noLowCoverage.49.conf %s/%s.%s.list.bz2 %s %sweightNbChr- %sweightNbChr+ -scoringMethod=%d" % (BaseRep,sign,ancF,anc.replace(' ','^'),weightM,weightP,scoring)
	
	if not os.access('%s.ok' % BaseName, os.R_OK):
		os.system('%s %s > %s.res.bz2 2> %s.log' % (Executable,Arguments,BaseName,BaseName))
		os.system('touch %s.ok' % BaseName)

