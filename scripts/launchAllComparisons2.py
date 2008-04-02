#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], "Lance la commande pour toutes les branches de l'arbre")

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.listSpecies):
	command = "/users/ldog/muffato/work/scripts/buildAncDiags.py %s -genesFile=data49/genes/genes.%%s.list.bz2 -ancGenesFile=data49/ancGenes/ancGenes.%%s.list.bz2 -target=%s,%s -showAncestral -showProjected 2>&1 | grep Extraction | sed 's/.*\.\.\. //' | sed 's/ OK//'" % (noms_fichiers["phylTree.conf"], e1.replace(' ', '^'), e2.replace(' ', '^'))
	print >> sys.stderr, e1, e2
	print "%s\t%s\t" % (e1,e2),
	sys.stdout.flush()
	os.system(command)
	sys.stdout.flush()


