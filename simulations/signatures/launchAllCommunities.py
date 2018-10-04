#!/usr/bin/env python2

import os
import sys
import itertools

import utils.myPhylTree

Executable = '/workspace/muffato/newSimu/getAncCommunitiesScores.py'

phylTree = utils.myPhylTree.PhylogeneticTree("/workspace/muffato/data49b/phylTree.noLowCoverage.49.conf")

sim = int(sys.argv[2])

for (anc,randomWalks,scoring,genes,weightM,weightP,threshold) in itertools.product(["Boreoeutheria"],[2,3,4,5,7,10,15,25],[0],"-+","-+","-+",[.01,.017,.02,.03]):


	ancF = phylTree.fileName[anc]
	BaseRep = "/workspace2/muffato/newSimu/all2/%d/" % sim
	BaseName = BaseRep + "scores/%s.X.%s.%d.%s.%s" % (ancF,genes,scoring,weightM,weightP)
	Arguments = "+psyco +bz2 /workspace/muffato/data49b/phylTree.noLowCoverage.49.conf %s/diags/%s.list -ancestr=%s +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ -scoringMethod=%d -walktrapLength=%d -genesFile=%s/genes/genes.%%s.list.bz2 -ancGenesFile=/workspace/muffato/data49b/ancGenes/ancGenes.%%s.list.bz2  -scoreThreshold=%f" % (BaseRep,ancF,anc.replace(' ','^'),genes,weightM,weightP,scoring,randomWalks,BaseRep, threshold)
	
	if not os.access('%s.ok' % BaseName, os.R_OK):
		os.system('%s %s > %s.res.bz2 2> %s.log' % (Executable,Arguments,BaseName,BaseName))
		os.system('touch %s.ok' % BaseName)

