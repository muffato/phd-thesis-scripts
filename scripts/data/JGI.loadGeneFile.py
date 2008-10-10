#! /users/ldog/muffato/python

__doc__ = """
	Telecharge depuis le site de l'UCSC le fichier avec les annotations xref et cree la liste des genes
"""

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("annotationFile",file))], \
	[("species",str,""), ("OUT.directory",str,""), \
	("OUT.genesFile",str,"genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"genes/full/genes.%s.list.bz2"), \
	("OUT.xrefFile",str,"xref/xref.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

OUTgenesFile = os.path.join(arguments["OUT.directory"], arguments["OUT.genesFile"])
OUTfullGenesFile = os.path.join(arguments["OUT.directory"], arguments["OUT.fullGenesFile"])
OUTxrefFile = os.path.join(arguments["OUT.directory"], arguments["OUT.xrefFile"])
for dir in [OUTgenesFile, OUTfullGenesFile, OUTxrefFile]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass

# On charge le fichier d'annotations
print >> sys.stderr, "Chargement du fichier et ecriture des genes ...",
foG = utils.myTools.myOpenFile(OUTgenesFile % phylTree.fileName[arguments["species"]], 'w')
foG2 = utils.myTools.myOpenFile(OUTfullGenesFile % phylTree.fileName[arguments["species"]], 'w')
dicXref = {}
for ligne in utils.myTools.myOpenFile(arguments["annotationFile"], 'r'):
	t = [intern(x) for x in ligne.replace('\n', '').split('\t')]
	print >> foG, "\t".join([t[x] for x in [2,4,5,3,1]])
	print >> foG2, "\t".join([t[x] for x in [2,4,5,3,1]])
	dicXref[t[1]] = t[-4]
foG.close()
foG2.close()
print >> sys.stderr, "%d genes OK" % len(dicXref)

# On va traiter les annotations maintenant

# Le dictionnaire xref
lstXref = set()
# D'abord charger celles deja connnues
print >> sys.stderr, "Chargement des annotations xref de reference ",
for esp in phylTree.listSpecies:
	for ligne in utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'r'):
		lstXref.update( ligne.replace('\n', '').split('\t')[3:] )
	sys.stderr.write('.')
print >> sys.stderr, " OK"

# On ecrit les notres seulement si elle sont deja connus
print >> sys.stderr, "Ecriture des xref ...",
foX = utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[arguments["species"]], 'w')
nb = 0
for (gene,xref) in dicXref.iteritems():
	reflink = list(lstXref.intersection(xref.split()))
	if len(reflink) > 0:
		nb += 1
		print >> foX, '\t'.join([gene,"-","-"]+reflink)
foX.close()
print >> sys.stderr, "%d xref OK" % nb

