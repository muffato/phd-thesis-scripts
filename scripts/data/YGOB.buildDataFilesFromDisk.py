#! /users/ldog/muffato/python


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import operator
import utils.myPhylTree
import utils.myTools


#############
# FONCTIONS #
#############

arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("inputDir",file), ("outputDir",file)], [], "Lit les donnees telechargees des levures et les formatte")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


dicStrand =  {"C": "-1", "W": "1"}
for esp in phylTree.listSpecies:
	continue
	f = open(arguments["inputDir"] + "/%s_genome.tab" % esp, "r")
	fout = utils.myTools.myOpenFile(arguments["outputDir"] + "/genes/genes.%s.list.bz2" % esp, "w")
	for (i,l) in enumerate(f):
		c = l.split()
		print >> fout, "\t".join( (c[1], str(i), str(i), dicStrand[c[2]], c[0]) )
	f.close()
	fout.close()


os.system("cat %s/Pillars.tab | sed 's/---//g' | sed 's/\t\+/ /g' | sed 's/^ //' | awk 'NF>1' | bzip2 > %s/ancGenes/ancGenes.Anc-Agos.list.bz2" % (arguments["inputDir"],arguments["outputDir"]))
os.system("cut --complement -f4 %s/Pillars.tab | sed 's/---//g' | sed 's/\t\+/ /g' | sed 's/^ //' | awk 'NF>1' | bzip2 > %s/ancGenes/ancGenes.Anc-Klac.list.bz2" % (arguments["inputDir"],arguments["outputDir"]))
os.system("cut --complement -f4,5 %s/Pillars.tab | sed 's/---//g' | sed 's/\t\+/ /g' | sed 's/^ //' | awk 'NF>1' | bzip2 > %s/ancGenes/ancGenes.Anc-PRE-DUP.list.bz2" % (arguments["inputDir"],arguments["outputDir"]))
os.system("cut -f6,7 %s/Pillars.tab | sed 's/---//g' | sed 's/\t\+/ /g' | sed 's/^ //' | awk 'NF>1' | bzip2 > %s/ancGenes/ancGenes.Anc-Skluy-Kwal.list.bz2" % (arguments["inputDir"],arguments["outputDir"]))


