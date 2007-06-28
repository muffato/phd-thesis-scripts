#! /users/ldog/muffato/python -OO

import os

#for seuilMaxDistInterGenes in [3,5,10,20,50]:
for seuilMaxDistInterGenes in [-1]:
	for useOutgroups in [0,1,2]:
		patt = "/users/ldog/muffato/workspace/tmp/simu.concorde.theria/%d.%d" % (seuilMaxDistInterGenes,useOutgroups)
		os.system( "/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildSortedGenome.py /users/ldog/muffato/workspace/tmp/simu.concorde/Genome.Theria.shuffled.bz2 /users/ldog/muffato/work/phylTree.vertebrates.44.conf -seuilMaxDistInterGenes=%d -useOutgroups=%d " % (seuilMaxDistInterGenes,useOutgroups) + "-ancestr=Theria -genesFile=/users/ldog/muffato/work/simu.constraint/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.constraint/ancGenes/ancGenes.%s.list.bz2 +withConcordeOutput / " + "%s.genome // %s.log"  % (patt,patt) )





