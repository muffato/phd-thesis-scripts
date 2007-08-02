#! /users/ldog/muffato/python -OO

import os

for seuilMaxDistInterGenes in [-1,3,5,10,20,50]:
#for seuilMaxDistInterGenes in [-1]:
	for useOutgroups in [0,1,2]:
		patt = "/users/ldog/muffato/workspace/tmp/concorde.theria.extreme/%d.%d" % (seuilMaxDistInterGenes,useOutgroups)
		os.system( "/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildSortedGenome.py /users/ldog/muffato/workspace/tmp/concorde.theria.extreme/Genome.Theria.shuffled.bz2 /users/ldog/muffato/work/phylTree.vertebrates.45.conf -seuilMaxDistInterGenes=%d -useOutgroups=%d " % (seuilMaxDistInterGenes,useOutgroups) + '-ancestr=Theria -genesFile=/users/ldog/muffato/work/simu.extreme/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.extreme/ancGenes/ancGenes.%s.list.bz2 +withConcordeOutput \\~ "(machine == \\\"dyogen04.ens.fr\\\")" / ' + "%s.genome // %s.log"  % (patt,patt) )





