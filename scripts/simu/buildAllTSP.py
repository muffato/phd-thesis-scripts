#!/usr/bin/env python2

import os
import itertools

print "Executable = /users/ldog/muffato/work/scripts/buildSortedGenome.py"
print "Universe = vanilla"
print "GetEnv = True"
print "Input = "
print "Initialdir = /users/ldog/muffato/workspace/tmp/allsimu"
print "should_transfer_files = NO"
print "Requirements = (target.Machine != \"heimdall.ens.fr\")"
print

for (seuilMaxDistInterGenes,useOutgroups) in itertools.product([-1,3,5,10,20,50],[0,1,2]):
	print "Output = %d/concorde/%d.%d.res" % (seuilMaxDistInterGenes,useOutgroups)
	print "Error = %d/concorde/%d.%d.log" % (seuilMaxDistInterGenes,useOutgroups)
	#print "Arguments = +psyco /users/ldog/muffato/work/phylTree.vertebrates.45.conf -fusionThreshold=%d -minimalLength=%d " % (fusionThreshold,minimalLength) + keepOnlyOrthos + "keepOnlyOrthos " + sameStrand + "sameStrand +useOutgroups -target=Euteleostomi +cutLongestPath -searchUndetectedSpecies -genesFile=%d/genes" % indSimu + "/genes.%s.list.bz2" + " -ancGenesFile=%d/ancGenes" % indSimu + "/ancGenes.%s.list.bz2"

	#print "/users/ldog/muffato/workspace/tmp/concorde.theria.extreme/Genome.Theria.shuffled.bz2 /users/ldog/muffato/work/phylTree.vertebrates.45.conf -seuilMaxDistInterGenes=%d -useOutgroups=%d " % (seuilMaxDistInterGenes,useOutgroups) + '-ancestr=Theria -genesFile=/users/ldog/muffato/work/simu.extreme/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.extreme/ancGenes/ancGenes.%s.list.bz2 +withConcordeOutput"

	#print "Log = /users/ldog/muffato/condor/alldiags.%d.%d.%d.%s.%s.condor_log" % (indSimu,fusionThreshold,minimalLength,keepOnlyOrthos,sameStrand)
	#print "Queue"
	#print







