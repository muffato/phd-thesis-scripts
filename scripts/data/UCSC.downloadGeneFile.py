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
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("xenoFile",file)], \
	[("species",str,""), ("speciesVersion",str,""), ("OUT.directory",str,""), \
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

# On telecharge le fichier de l'UCSC
combGenes = utils.myTools.defaultdict(list)
print >> sys.stderr, "Telechargement du fichier ...",
for ligne in utils.myTools.myOpenFile(arguments["xenoFile"], 'r'):
	t = ligne.replace('\n', '').split('\t')
	combGenes[t[0]].append( (int(t[1]),int(t[2]), t[3], set(t[4:])) )

# Le dictionnaire xref
lstXref = set()
# D'abord charger celles deja connnues
print >> sys.stderr, "Chargement des annotations xref de reference ",
for esp in phylTree.listSpecies:
	if phylTree.officialName[arguments["species"]] == esp:
		continue
	if not utils.myTools.fileAccess(OUTxrefFile % phylTree.fileName[esp]):
		continue
	for ligne in utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'r'):
		lstXref.update( ligne.replace('\n', '').split('\t')[3:] )
	sys.stderr.write('.')
print >> sys.stderr, " OK"


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
		# On ecrit les notres seulement si elle sont deja connus
		genes.append( (c,x1,x2,strand,list(reflink.intersection(lstXref))) )
		#genes.append( (c,x1,x2,strand,reflink) )
del combGenes
noms = ["MYGENE.%s.%06d" % (arguments["speciesVersion"],i+1) for i in xrange(len(genes))]
print >> sys.stderr, "%d genes OK" % len(genes)

esp = phylTree.fileName[arguments["species"]]

# Ecriture des fichiers de genes / xref
print >> sys.stderr, "Ecriture des fichiers ...",

foG = utils.myTools.myOpenFile(OUTgenesFile % esp, 'w')
foX = utils.myTools.myOpenFile(OUTxrefFile % esp, 'w')
for (i,(c,x1,x2,strand,reflink)) in enumerate(genes):
	print >> foG, "\t".join([str(x) for x in [c,x1,x2,strand,noms[i]]])
	print >> foX, '\t'.join([noms[i],"-","-"]+list(reflink))
foG.close()
foX.close()

os.symlink(OUTgenesFile % esp, OUTfullGenesFile % esp)

print >> sys.stderr, "OK"

