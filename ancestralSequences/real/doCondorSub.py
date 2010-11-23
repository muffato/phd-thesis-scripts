#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("lstCDNA",file), ("treesFile",file)], [("genesFile",str,"")], "Cree le fichier de soumission a condor")

# Association "AAEL000042-RA" > "./Aedes_aegypti/AAEL000042-RA.fa"
dic = {}
f = utils.myFile.openFile(arguments["lstCDNA"], "r")
for l in f:
	dic[ l[:-4].split('/')[2] ] = l[:-1]
f.close()

# Association gene_name > transcript_name
phylTree = utils.myPhylTree.PhylogeneicTree(arguments["phylTree.conf"])
notransc = {}
for e in phylTree.listSpecies:
	f = utils.myFile.openFile(arguments["genesFile"] % phylTree.fileName[e], "r")
	for l in f:
		t = l.split("\t")
		if len(t) >= 7:
			notransc[t[4]] = t[6]
	f.close()

print """
Executable = /usr/bin/python2.5
Universe = vanilla
Input =
Log = /workspace/muffato/condorlog
GetEnv = True
should_transfer_files = NO
Requirements = (target.Machine != "heimdall.ens.fr")
Notify_user = muffato@biologie.ens.fr
Notification = complete
NiceUser = True
Rank = kflops
"""


f = utils.myFile.openFile(arguments["treesFile"], "r")
for (i,l) in enumerate(f):
	names = l.replace("(", " ").replace(")", " ").replace(",", " ").split()
	# Le cas de la vache dans Ensembl53
	names = [s.split("/")[0] for s in names]
	names = [s if s in dic else notransc[s] for s in names]

	print "Initialdir = /users/ldog/muffato/workspace/ancestralSequences/real/fam/%d" % (i+1)
	print "Output = RES"
	print "Error = LOG"
	print "Arguments =  ../../../ortheus/Ortheus.py -j -d %s -e %s" % (l[:-1], " ".join("seq/%s.fa" % x for x in names))
	print "Queue"
	print
	continue

	for s in names:
		s = s.split("/")[0]
		if s not in dic:
			s = notransc[s]
			#print >> sys.stderr, s
			#continue
		print "ln -s ../../../cds/%s fam/%d/seq/%s" % (dic[s], i+1, s)
f.close()


