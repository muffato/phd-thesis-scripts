#! /users/ldog/muffato/python -OO

__doc__ = """
	Calcule les taux de GC des especes modernes
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("CDSfile",str,""), ("GCfile",str,"")], __doc__ )


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Le repertoire
try:
	os.makedirs(os.path.dirname(options["GCfile"]))
except OSError:
	pass

for esp in phylTree.listSpecies:
	
	print >> sys.stderr, "Parsing %s ..." % esp,
	fi = utils.myTools.myOpenFile(options["CDSfile"] % phylTree.fileName[esp], 'r')
	fo = utils.myTools.myOpenFile(options["GCfile"] % phylTree.fileName[esp], 'w')
	for l in fi:
		(gene,seq) = l[:-1].split("\t")
		seq = seq.upper()
		dicGC = [0, 0, 0]
		dicALL = [0, 0, 0]
		for (i,base) in enumerate(seq):
			dicALL[i % 3] += 1
			if base in "GC":
				dicGC[i % 3] += 1
		res = [gene] + dicGC + dicALL + [(100.*dicGC[i])/dicALL[i] for i in xrange(3)] + [100.*float(sum(dicGC)-dicGC[2])/float(sum(dicALL)-dicALL[2]), float(100*sum(dicGC))/sum(dicALL)]
		print >> fo, utils.myTools.printLine(res)
	fo.close()
	fi.close()
	print >> sys.stderr, "OK"

