#! /users/ldog/muffato/python -OO

import os
import utils.myTools

print "Executable = /users/ldog/muffato/work/scripts/buildAncCommunities.py"
print "Universe = vanilla"
print "GetEnv = True"
print "Input = "
print "Initialdir = /users/ldog/muffato/workspace/tmp/communities.boreo.extreme/"
print "should_transfer_files = NO"
print "Requirements = (target.Machine != \"heimdall.ens.fr\")"


for (weightM,weightP,genes,distScoring,randomWalks) in utils.myTools.myIterator.tupleOnManyLists("+-","+-","+-","+-",[2,3,4,5,7,10,15,25]):
	print "Output = /users/ldog/muffato/workspace/tmp/communities.boreo.extreme/%d.%s.%s.%s.%s.res" % (randomWalks,distScoring,genes,weightM,weightP)
	print "Error = /users/ldog/muffato/workspace/tmp/communities.boreo.extreme/%d.%s.%s.%s.%s.log" % (randomWalks,distScoring,genes,weightM,weightP)
	print "Arguments = +psyco /users/ldog/muffato/work/phylTree.vertebrates.45.conf /workspace/muffato/tmp/diags.all.extreme/0.2.-.+.all.res.bz2 -ancestr=Boreoeutheria +useOutgroups " + genes + "useLonelyGenes " + weightM + "weightNbChr- " + weightP + "weightNbChr+ " + distScoring + "newScoring -walktrapLength=%d" % (randomWalks) + " -genesFile=/users/ldog/muffato/work/simu.extreme/genes/genes.%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/simu.extreme/ancGenes/ancGenes.%s.list.bz2"
	print "Log = /users/ldog/muffato/condor/communities.boreo.extreme.%d.%s.%s.%s.%s.condor_log" % (randomWalks,distScoring,genes,weightM,weightP)
	print "Queue"

