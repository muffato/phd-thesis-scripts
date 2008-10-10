#! /users/ldog/muffato/python

import os
import itertools

print "Executable = /users/ldog/muffato/work/scripts/buildAncDiags.py"
print "Universe = vanilla"
print "GetEnv = True"
print "Input = "
print "Initialdir = /users/ldog/muffato/workspace/tmp/allsimu"
print "should_transfer_files = NO"
print "Requirements = (target.Machine != \"heimdall.ens.fr\")"
print

for (indSimu,keepOnlyOrthos,fusionThreshold,minimalLength,sameStrand) in itertools.product([1,2,3,4,5],"+-",[-1,0,1,2],[2,3,4],"+-"):
	print "Output = %d/diags/%d.%d.%s.%s.res" % (indSimu,fusionThreshold,minimalLength,keepOnlyOrthos,sameStrand)
	print "Error = %d/diags/%d.%d.%s.%s.log" % (indSimu,fusionThreshold,minimalLength,keepOnlyOrthos,sameStrand)
	print "Arguments = +psyco /users/ldog/muffato/work/phylTree.vertebrates.45.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength) + keepOnlyOrthos + "keepOnlyOrthos " + sameStrand + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=%d/genes" % indSimu + "/genes.%s.list.bz2" + " -ancGenesFile=%d/ancGenes" % indSimu + "/ancGenes.%s.list.bz2"


	print "Log = /users/ldog/muffato/condor/alldiags.%d.%d.%d.%s.%s.condor_log" % (indSimu,fusionThreshold,minimalLength,keepOnlyOrthos,sameStrand)
	print "Queue"
	print

