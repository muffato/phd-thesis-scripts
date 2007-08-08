#! /users/ldog/muffato/python -OO

__doc__ = """
	Cree des fichiers de donnees vide
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
	("OUT.genesFile",str,"genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"genes/full/genes.%s.list.bz2"), \
	("OUT.xrefFile",str,"xref/xref.%s.list.bz2"), \
	("OUT.orthosFile",str,"orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.paras2File",str,"orthologs/paras.%s.%s.list.bz2"), \
	("OUT.paras1File",str,"paralogs/paras.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTgenesFile = os.path.join(options["OUT.directory"], options["OUT.genesFile"])
OUTfullGenesFile = os.path.join(options["OUT.directory"], options["OUT.fullGenesFile"])
OUTxrefFile = os.path.join(options["OUT.directory"], options["OUT.xrefFile"])
OUTorthosFile = os.path.join(options["OUT.directory"], options["OUT.orthosFile"])
OUTparas2File = os.path.join(options["OUT.directory"], options["OUT.paras2File"])
OUTparas1File = os.path.join(options["OUT.directory"], options["OUT.paras1File"])
for dir in [OUTgenesFile, OUTfullGenesFile, OUTxrefFile, OUTorthosFile, OUTparas1File, OUTparas2File]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass


def touch(nom):
	if not os.access(nom, os.F_OK):
		print >> sys.stderr, nom
		f = utils.myTools.myOpenFile(nom, 'w')
		f.close()

for esp in phylTree.listSpecies:
	
	touch(OUTgenesFile % phylTree.fileName[esp])
	touch(OUTfullGenesFile % phylTree.fileName[esp])
	touch(OUTxrefFile % phylTree.fileName[esp])
	touch(OUTparas1File % phylTree.fileName[esp])

# On genere les fichiers d'homologues
for (esp1,esp2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.listSpecies):

	touch(OUTorthosFile % (phylTree.fileName[esp1],phylTree.fileName[esp2]))
	touch(OUTorthosFile % (phylTree.fileName[esp2],phylTree.fileName[esp1]))
	touch(OUTparas2File % (phylTree.fileName[esp1],phylTree.fileName[esp2]))
	touch(OUTparas2File % (phylTree.fileName[esp2],phylTree.fileName[esp1]))

