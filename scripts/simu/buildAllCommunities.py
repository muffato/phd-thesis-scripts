#! /users/ldog/muffato/python

import os
import sys
import itertools
import utils.myPhylTree

Executable = '/users/ldog/muffato/work/scripts/buildAncCommunities.py'

phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

for (anc,randomWalks,scoring,genes,weightM,weightP) in itertools.product([sys.argv[2]],[2,3,4,5,7,10,15,25],[0,1,2],"-+","-+","-+"):


	ancF = phylTree.fileName[anc]
	BaseName = "/workspace/muffato/ancGenomes47/%s.%d.%s.%d.%s.%s" % (ancF,randomWalks,genes,scoring,weightM,weightP)
	Arguments = "/users/ldog/muffato/work/data47/phylTree.vertebrates.47.conf /users/ldog/muffato/work/data47/diags/anc/%s.list.bz2 -ancestr=%s +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ -scoringMethod=%d -walktrapLength=%d -genesFile=/users/ldog/muffato/work/data47/genes/genes.%%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/data47//ancGenes/ancGenes.%%s.list.bz2" % (ancF,anc.replace(' ','^'),genes,weightM,weightP,scoring,randomWalks)
	
	if not os.access('%s.ok' % BaseName, os.R_OK):
		os.system('%s %s > %s.res 2> %s.log' % (Executable,Arguments,BaseName,BaseName))
		os.system('touch %s.ok' % BaseName)

