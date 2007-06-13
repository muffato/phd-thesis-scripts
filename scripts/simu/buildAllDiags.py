#! /users/ldog/muffato/python -OO

import os

transform = {True:"+", False:"-"}

for keepOnlyOrthos in [True,False]:
	for fusionThreshold in [-1,0,1,2]:
		for minimalLength in [2,3,4]:
			for sameStrand in [True,False]:
				patt = "/users/ldog/muffato/workspace/%d.%d." % (fusionThreshold,minimalLength) + transform[keepOnlyOrthos] + transform[sameStrand]
				os.system( ("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncDiags.py /users/ldog/muffato/work/phylTree.vertebrates.44.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength)) + transform[keepOnlyOrthos] + "keepOnlyOrthos " + transform[sameStrand] + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=/users/ldog/muffato/work/simu.constraint/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.constraint/ancGenes/ancGenes.%s.list.bz2 / " + "%s.all.res // %s.all.log"  % (patt,patt) )
				os.system( ("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncDiags.py /users/ldog/muffato/work/phylTree.noLowCoverage.44.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength)) + transform[keepOnlyOrthos] + "keepOnlyOrthos " + transform[sameStrand] + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=/users/ldog/muffato/work/simu.constraint/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.constraint/ancGenes/ancGenes.%s.list.bz2 / " + "%s.noLow.res // %s.noLow.log" % (patt,patt) )





