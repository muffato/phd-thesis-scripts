#! /users/ldog/muffato/python

__doc__ = """
	Cree le tableau des stats des blocs de syntenie et des genes ancestraux en fonction d'un seuil de coupure
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("ancGenesFile",str,"ancGenes/ancGenes.%s.list.bz2"), ("diagsFile",str,"diags/anc/diags.%s.list.bz2"), ("logFile",str,"log.diags")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Liste des especes dans le bon ordre
todo = set(phylTree.listAncestr)
l1 = phylTree.dicLinks["Euteleostomi"]["Homo sapiens"][:-1]
todo.difference_update(l1)
l2 = phylTree.dicLinks["Glires"]["Murinae"]
todo.difference_update(l2)
l3 = [e for e in todo if phylTree.isChildOf(e, "Mammalia")]
l3 = sorted(l3, key = lambda e: phylTree.ages[e], reverse = True)
todo.difference_update(l3)
l4 = [e for e in todo if phylTree.isChildOf(e, "Clupeocephala")]
l4 = sorted(l4, key = lambda e: phylTree.ages[e], reverse = True)
todo.difference_update(l4)
l5 = [e for e in todo if phylTree.isChildOf(e, "Amniota")]
l5 = sorted(l5, key = lambda e: phylTree.ages[e], reverse = True)
todo.difference_update(l5)
l6 = sorted(todo, key = lambda e: phylTree.ages[e], reverse = True)
lstEspeces = l6 + l5 + l4 + l1 + l3 + l2
#lstEspeces = l5

ref = None

for cutoff in ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"]:

	# Nom et age
	data = {}
	for e in lstEspeces:
		data[e] = [e, phylTree.ages[e]]

	# Nombre de genes
	for e in lstEspeces:
		data[e].append( len(utils.myGenomes.Genome(cutoff + "/" + (arguments["ancGenesFile"] % phylTree.fileName[e])).lstGenes[None]) )

	# Nombre de blocs
	for e in lstEspeces:
		f = utils.myFile.openFile(cutoff + "/" + (arguments["diagsFile"] % phylTree.fileName[e]), "r")
		lst = []
		for l in f:
			x = int(l.split("\t")[1])
			if x >= 2:
				lst.append(x)
		f.close()
		data[e].append(len(lst))
		data[e].append(sum(lst))
		data[e].append((100.*data[e][4])/data[e][2])

	# Stats sur les blocs
	f = utils.myFile.openFile(cutoff + "/" + arguments["logFile"], "r")
	for l in f:
		if "ancestral" not in l:
			continue
		i = l.index(".")
		e = l[41:i-1]
		s = l[i+4:-4].replace("/"," ").replace("["," ").replace("]"," ").replace("-"," ").split()[:11]
		if e in data:
			data[e].extend(s)
	f.close()

	print utils.myFile.myTSV.printLine(["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean", "Stddev"])
	for e in lstEspeces:
		print utils.myFile.myTSV.printLine(data[e])

	if ref is None:
		ref = data
	else:
		print utils.myFile.myTSV.printLine(["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "%Useful Gene Loss", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean"])
		def trans(x):
			if type(x) == str:
				try:
					y = int(x)
					return y
				except:
					try:
						y = float(x)
						return y
					except:
						return x
			else:
				return x
		for e in lstEspeces:
			newdata = [(trans(x)-trans(ref[e][i]) if i >= 2 else x) for (i,x) in enumerate(data[e][:-1])]
			newdata.insert(6, 100*(1.-float(newdata[4])/newdata[2]) if newdata[2] != 0 else None)
			print utils.myFile.myTSV.printLine(newdata)


