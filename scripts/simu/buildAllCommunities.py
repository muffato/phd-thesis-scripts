#! /users/ldog/muffato/python

import os
import sys
import itertools
import utils.myPhylTree

Executable = '/users/ldog/muffato/work/scripts/buildAncCommunities2.py'


for (anc,randomWalks,scoring,genes,weightM,weightP) in itertools.product([sys.argv[2]],[2,3,4,5,7,10,15,25],[0,1,2,3,4],[1,2],"-+","-+"):


	ancF = phylTree.fileName[anc]
	BaseName = "/workspace/muffato/ancGenomes47/%s.%d.%s.%d.%s.%s" % (ancF,randomWalks,genes,scoring,weightM,weightP)
	Arguments = "/users/ldog/muffato/work/data47/phylTree.vertebrates.47.conf /users/ldog/muffato/work/data47/diags/anc/%s.list.bz2 -ancestr=%s +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ -scoringMethod=%d -walktrapLength=%d -genesFile=/users/ldog/muffato/work/data47/genes/genes.%%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/data47//ancGenes/ancGenes.%%s.list.bz2" % (ancF,anc.replace(' ','^'),genes,weightM,weightP,scoring,randomWalks)
	
	if not os.access('%s.ok' % BaseName, os.R_OK):
		os.system('%s %s > %s.res 2> %s.log' % (Executable,Arguments,BaseName,BaseName))
		os.system('touch %s.ok' % BaseName)


for (randomWalks,scoring,genes,weightM,weightP) in itertools.product([2,3,4,5,7,10,15,25],range(5),[1,2],"-+","-+"):

	BaseName = "scores/%s/raw/Boreo.%d.%d.%s.%s" % (sys.argv[1], genes,scoring,weightM,weightP)
	Arguments = "../data51/phylTree.%s.51.conf diagsHuman Boreoeutheria -genesFile=../data51/genes/genes.%%s.list.bz2 -ancGenesFile=../data51/ancGenes/ancGenes.%%s.list.bz2 +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ -scoringMethod=%d" % (sys.argv[2], genes,weightM,weightP,scoring)

	os.system('condor-submit.sh ./buildAncCommunities.py %s / %s.res // %s.log' % (Arguments,BaseName,BaseName))

for (randomWalks,scoring,genes,weightM,weightP) in itertools.product([2,3,4,5,7,10,15,25],range(5),[1,2],"-+","-+"):
#for (randomWalks,scoring,genes,weightM,weightP) in itertools.product([2],[0],[1],"-","-"):

	BaseName = "genomes/%s/Boreo.%d.%d.%s.%s" % (sys.argv[1],genes,scoring,weightM,weightP)
	Arguments = "../data51/phylTree.%s.51.conf ../data51/diags/full/anc/diags.Boreoeutheria.list.bz2 -minDiagLength=%d %sweightNbChr- %sweightNbChr+ -scoringMethod=%d -weightFile=../data51/weights.lst" % (sys.argv[1], genes,weightM,weightP,scoring)
	if scoring == 2:
		Arguments = Arguments + " +withMemoize"

	os.system('condor-submit.sh ../scripts/buildAncCommunities2.py %s / %s.res // %s.log $' % (Arguments,BaseName,BaseName))

