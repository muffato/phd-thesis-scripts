#! /users/ldog/muffato/python -OO

import os


for weightM in "+-":
	for weightP in "+-":
		for genes in "+-":
			for distScoring in "+-":
				for randomWalks in [2,3,4,5,7,10,15,25]:
					for qualityFunction in [1,2,3]:
						patt = "/users/ldog/muffato/workspace/%d.%d.%s.%s.%s.%s" % (qualityFunction,randomWalks,distScoring,genes,weightM,weightP)
						os.system("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncCommunities.py /users/ldog/muffato/work/phylTree.vertebrates.44.conf /workspace/muffato/tmp/simu.walktrap/res.diags -ancestr=Boreoeutheria +useOutgroups " + genes + "useLonelyGenes " + weightM + "weightNbChr- " + weightP + "weightNbChr+ " + distScoring + "newScoring -walktrapLength=%d -qualityFunction=%d" % (randomWalks,qualityFunction) + " -genesFile=/users/ldog/muffato/work/simu.constraint2/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.constraint2/ancGenes/ancGenes.%s.list.bz2 / " + "%s.all.res // %s.all.log"  % (patt,patt) )
						




