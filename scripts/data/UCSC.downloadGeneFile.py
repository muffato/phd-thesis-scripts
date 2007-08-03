#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site de l'UCSC le fichier avec les annotations xref et cree la liste des genes
"""

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("species",str,""), ("database",str,""), ("OUT.directory",str,""), \
	("IN.UCSC-URL",str,"http://hgdownload.cse.ucsc.edu/goldenPath/XXX/database/xenoRefFlat.txt.gz"), \
	("OUT.genesFile",str,"genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"full/genes.%s.list.bz2"), \
	("OUT.xrefFile",str,"xref/xref.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTgenesFile = os.path.join(options["OUT.directory"], options["OUT.genesFile"])
OUTfullGenesFile = os.path.join(options["OUT.directory"], options["OUT.fullGenesFile"])
OUTxrefFile = os.path.join(options["OUT.directory"], options["OUT.xrefFile"])
for dir in [OUTgenesFile, OUTfullGenesFile, OUTxrefFile]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass

# On telecharge le fichier de l'UCSC
combGenes = utils.myTools.defaultdict(list)
print >> sys.stderr, "Telechargement du fichier ...",
for ligne in utils.myTools.myOpenFile(options["IN.UCSC-URL"].replace("XXX", options["database"]), 'r'):
	t = [intern(x) for x in ligne[:-1].split('\t')]
	combGenes[t[2]].append( (int(t[4]),int(t[5]),t[3], set([t[0],t[1]]) ) )

# On prend les xref chevauchants
print >> sys.stderr, "Combinaison des alignements ...",
genes = []
for (c,comb) in combGenes.iteritems():
	comb.sort()
	while len(comb) > 0:
		(x1,x2,strand,reflink) = comb.pop(0)
		while len(comb) > 0:
			(x1B,x2B,_,reflinkB) = comb[0]
			if x1B <= x2:
				x2 = max(x2, x2B)
				reflink.update(reflinkB)
				del comb[0]
			else:
				break
		genes.append( (c,x1,x2,strand,list(reflink)) )
del combGenes
noms = ["MYGENE.%s.%06d" % (options["database"],i+1) for i in xrange(len(genes))]
print >> sys.stderr, "%d genes OK" % len(genes)

esp = options["species"]

# Ecriture des fichiers de genes / xref
print >> sys.stderr, "Ecriture des fichiers ...",

foG = utils.myTools.myOpenFile(OUTgenesFile % phylTree.fileName[esp], 'w')
foG2 = utils.myTools.myOpenFile(OUTgenesFile % phylTree.fileName[esp], 'w')
foX = utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'w')
for (i,(c,x1,x2,strand,reflink)) in enumerate(genes):
	print >> foG, "\t".join([str(x) for x in [c,x1,x2,strand,noms[i]]])
	print >> foG2, "\t".join([str(x) for x in [c,x1,x2,strand,noms[i]]])
	print >> foX, '\t'.join([noms[i],"-","-"]+reflink)
foG.close()
foG2.close()
foX.close()

print >> sys.stderr, "OK"

