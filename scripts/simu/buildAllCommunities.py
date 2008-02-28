#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools
import utils.myPhylTree

print 'Executable = /users/ldog/muffato/work/scripts/buildAncCommunities.py'
print 'Universe = vanilla'
print 'GetEnv = True'
print 'Input = '
print 'Initialdir = /workspace/muffato/ancGenomes47'
print 'should_transfer_files = NO'
print 'Requirements = (target.Machine != "heimdall.ens.fr") && (target.Machine != "dyogen04.ens.fr")'
print 'Notify_user = muffato@biologie.ens.fr'
print


phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])
#indSim = int(sys.argv[2])

#for (anc,randomWalks,scoring,genes,weightM,weightP) in utils.myTools.myIterator.tupleOnManyLists(phylTree.listAncestr,[2,3,4,5,7,10,15,25],[0,1,2],"-+","-+","-+"):
for (anc,randomWalks,scoring,genes,weightM,weightP) in utils.myTools.myIterator.tupleOnManyLists([sys.argv[2]],[2,3,4,5,7,10,15,25],[0,1,2],"-+","-+","-+"):
	ancF = phylTree.fileName[anc]
	print "Output = /workspace/muffato/ancGenomes47/%s.%d.%s.%d.%s.%s.res" % (ancF,randomWalks,genes,scoring,weightM,weightP)
	print "Error = /workspace/muffato/ancGenomes47/%s.%d.%s.%d.%s.%s.log" % (ancF,randomWalks,genes,scoring,weightM,weightP)
	print "Arguments = /users/ldog/muffato/work/data47/phylTree.vertebrates.47.conf /users/ldog/muffato/work/data47/diags/anc/%s.list.bz2 -ancestr=%s +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ -scoringMethod=%d -walktrapLength=%d -genesFile=/users/ldog/muffato/work/data47/genes/genes.%%s.list.bz2 -ancGenesFile=/users/ldog/muffato/work/data47//ancGenes/ancGenes.%%s.list.bz2" % (ancF,anc.replace(' ','^'),genes,weightM,weightP,scoring,randomWalks)
	print "Log = /users/ldog/muffato/condor/communities.%s.%d.%s.%d.%s.%s.condor_log" % (ancF,randomWalks,genes,scoring,weightM,weightP)
	print "Queue"
	print

for (anc,randomWalks,distScoring,genes,weightM,weightP) in utils.myTools.myIterator.tupleOnManyLists(phylTree.listAncestr,[2,3,4,5,7,10,15,25],"-+","-+","-+","-+"):
	continue
	ancF = phylTree.fileName[anc]
	print "Output = /workspace/muffato/tmp/allsimu/%d/communities/%s.%d.%s.%s.%s.%s.res" % (indSim,ancF,randomWalks,distScoring,genes,weightM,weightP)
	print "Error = /workspace/muffato/tmp/allsimu/%d/communities/%s.%d.%s.%s.%s.%s.log" % (indSim,ancF,randomWalks,distScoring,genes,weightM,weightP)
	print "Arguments = +psyco /users/ldog/muffato/work/phylTree.vertebrates.45.conf /workspace/muffato/tmp/allsimu/%d/res.diags.bz2 -ancestr=%s +useOutgroups %suseLonelyGenes %sweightNbChr- %sweightNbChr+ %snewScoring -walktrapLength=%d -genesFile=/workspace/muffato/tmp/allsimu/%d/genes/genes.%%s.list.bz2 -ancGenesFile=/workspace/muffato/tmp/allsimu/%d/ancGenes/ancGenes.%%s.list.bz2" % (indSim,anc.replace(' ','^'),genes,weightM,weightP,distScoring,randomWalks,indSim,indSim)
	print "Log = /users/ldog/muffato/condor/communities.%d.%s.%d.%s.%s.%s.%s.condor_log" % (indSim,ancF,randomWalks,distScoring,genes,weightM,weightP)
	print "Queue"
	print

