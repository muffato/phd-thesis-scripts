#! /users/ldog/muffato/python -OO

__doc__ = """
	Calcule les taux de GC des especes ancestrales
"""

import sys
import utils.myTools
import utils.myPhylTree

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("ancGenesFile",str,""), ("GCmFile",str,""), ("GCaFile",str,"")], __doc__ )

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

utils.myTools.mkDir(options["GCaFile"])

dicALLGC = {}
dicSpecies = {}

for esp in phylTree.listSpecies:
	
	print >> sys.stderr, "Parsing %s ..." % esp,
	fi = utils.myTools.myOpenFile(options["GCmFile"] % phylTree.fileName[esp], 'r')
	for l in fi:
		t = l[:-1].split("\t")
		dicALLGC[t[0]] = tuple([int(x) for x in t[1:7]] + [float(x) for x in t[7:]])
		dicSpecies[t[0]] = esp
	fi.close()
	print >> sys.stderr, "OK"

for anc in phylTree.listAncestr:
	print >> sys.stderr, "Calculating ancestral %s ..." % anc,
	fi = utils.myTools.myOpenFile(options["ancGenesFile"] % phylTree.fileName[anc], 'r')
	fo = utils.myTools.myOpenFile(options["GCaFile"] % phylTree.fileName[anc], 'w')
	phylTree.initCalcDist(anc, False)
	for l in fi:
		t = l.split()
	
		# Les simples sommes et moyennes
		res = [sum([dicALLGC[x][i] for x in t]) for i in xrange(6)]
		if 0 in res[3:]:
			res += [0,0,0,0,0]
		else:
			res += [(100.*res[i])/res[i+3] for i in xrange(3)] + [100.*float(res[0]+res[1])/(res[3]+res[4]),100.*float(res[0]+res[1]+res[2])/(res[3]+res[4]+res[5])]
		
		# Les valeurs selon l'arbre
		for i in xrange(6,11):
			valL = utils.myTools.defaultdict(list)
			for x in t:
				valL[dicSpecies[x]].append(dicALLGC[x][i])
			values = {}
			for e in valL:
				values[e] = sum(valL[e]) / len(valL[e])
			res.append(phylTree.calcDist(values))
		print >> fo, utils.myTools.printLine(res)
	fo.close()
	fi.close()
	print >> sys.stderr, "OK"


