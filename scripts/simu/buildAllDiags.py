#! /users/ldog/muffato/python -OO

import os


for keepOnlyOrthos in "+-":
	for fusionThreshold in [-1,0,1,2]:
		for minimalLength in [2,3,4]:
			for sameStrand in "+-":
				patt = "/workspace/muffato/tmp/diags.simu.dup.%d.%d.%s.%s" % (fusionThreshold,minimalLength,keepOnlyOrthos,sameStrand)
				os.system( ("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncDiags.py /users/ldog/muffato/work/phylTree.vertebrates.45.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength)) + keepOnlyOrthos + "keepOnlyOrthos " + sameStrand + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=/workspace/muffato/tmp/genes.simu.dup/genes/genes.%s.list.bz2 -ancGenesFile=/workspace/muffato/tmp/genes.simu.dup/ancGenes/ancGenes.%s.list.bz2 / " + "%s.all.res // %s.all.log"  % (patt,patt) )
				os.system( ("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncDiags.py /users/ldog/muffato/work/phylTree.noLowCoverage.45.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength)) + keepOnlyOrthos + "keepOnlyOrthos " + sameStrand + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=/workspace/muffato/tmp/genes.simu.dup/genes/genes.%s.list.bz2 -ancGenesFile=/workspace/muffato/tmp/genes.simu.dup/ancGenes/ancGenes.%s.list.bz2 / " + "%s.noLow.res // %s.noLow.log"  % (patt,patt) )
				#os.system( ("/users/ldog/muffato/work/scripts/condor-submit.sh /users/ldog/muffato/work/scripts/buildAncDiags.py /users/ldog/muffato/work/phylTree.noLowCoverage.44.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength)) + keepOnlyOrthos + "keepOnlyOrthos " + sameStrand + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=/users/ldog/muffato/work/simu.constraint/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.constraint/ancGenes/ancGenes.%s.list.bz2 / " + "%s.noLow.res // %s.noLow.log" % (patt,patt) )





