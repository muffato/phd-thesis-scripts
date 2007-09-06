#! /users/ldog/muffato/python -OO

__doc__ = """
	Cree des liens entre orthos.A.B et orthos.B.A (idem paralogues)
"""

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("OUT.directory",str,""), \
	("OUT.orthosFile",str,"orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.paras2File",str,"orthologs/paras.%s.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTorthosFile = os.path.join(options["OUT.directory"], options["OUT.orthosFile"])
OUTparas2File = os.path.join(options["OUT.directory"], options["OUT.paras2File"])

# On cree les liens symboliques
for (esp1,esp2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.listSpecies):

	os.unlink(OUTorthosFile % (phylTree.fileName[esp2],phylTree.fileName[esp1]))
	os.symlink(os.path.basename(OUTorthosFile % (phylTree.fileName[esp1],phylTree.fileName[esp2])), OUTorthosFile % (phylTree.fileName[esp2],phylTree.fileName[esp1]))

	os.unlink(OUTparas2File % (phylTree.fileName[esp2],phylTree.fileName[esp1]))
	os.symlink(os.path.basename(OUTparas2File % (phylTree.fileName[esp1],phylTree.fileName[esp2])), OUTparas2File % (phylTree.fileName[esp2],phylTree.fileName[esp1]))

