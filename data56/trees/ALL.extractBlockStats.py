#!/usr/bin/env python2

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

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("diagsFile",str,"diags/anc/diags.%s.list.bz2"), ("outputODS",str,"")], \
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

#lstEspeces = ["Euteleostomi", "Amniota", "Boreoeutheria"]

allCutoff = ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50", "0.xx"]
#allCutoff = ["0.00", "0.05", "0.10"]
allCutoff = ["v1/x.04"]

alldata = {}
alldiff = {}

for cutoff in allCutoff:

	print >> sys.stderr, cutoff, "...",

	# Recuperation des donnees de longueur de blocs
	alldata[cutoff] = data = {}
	for e in lstEspeces:

		f = utils.myFile.openFile(cutoff + "/" + (arguments["diagsFile"] % phylTree.fileName[e]), "r")
		lst = []
		sing = 0
		tot = 0
		interv = 0
		for l in f:
			x = int(l.split("\t")[1])
			tot += x
			if x >= 2:
				lst.append(x)
				interv += (x-1)
			else:
				sing += 1
		f.close()

		data[e] = [e, phylTree.ages[e], tot, len(lst), tot-sing, (100.*(tot-sing))/tot, interv, (100.*interv)/(tot-20.)]
		data[e].extend(utils.myMaths.myStats.valSummary(lst)[:-1])
	
	if cutoff == allCutoff[0]:
		ref = data
	else:

		alldiff[cutoff] = diff = {}
		for e in lstEspeces:
			newdata = [(x-ref[e][i] if i >= 2 else x) for (i,x) in enumerate(data[e][:-1])]
			newdata.insert(6, 100*(1.-float(newdata[4])/newdata[2]) if newdata[2] != 0 else None)
			diff[e] = newdata
	print >> sys.stderr, "OK"

if arguments["outputODS"] == "":
	for cutoff in allCutoff:
		print utils.myFile.myTSV.printLine(["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean", "Stddev"])
		for e in lstEspeces:
			print utils.myFile.myTSV.printLine(alldata[cutoff][e])
	if cutoff in alldiff:
		print utils.myFile.myTSV.printLine(["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "%Useful Gene Loss", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean"])
		for e in lstEspeces:
			print utils.myFile.myTSV.printLine(alldiff[cutoff][e])

else:
	import odf.opendocument
	import odfpy.datatable

	textdoc = odf.opendocument.OpenDocumentSpreadsheet()

	for cutoff in allCutoff:

		# Premiere table avec les stats brutes
		val = [["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "NbInt", "%CovInt", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean", "Stddev"]]
		for e in lstEspeces:
			val.append(alldata[cutoff][e])

		table = odfpy.datatable.DataTable(val)
		table.datasourcehaslabels = "both"
		t = table()
		t.setAttribute('name', cutoff)
		textdoc.spreadsheet.addElement(t)

		if cutoff in alldiff:

			# Deuxieme table avec les differences par rapport a la reference
			val = [["Ancestor", "Age (My)", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "%Useful Gene Loss", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean"]]
			for e in lstEspeces:
				val.append(alldiff[cutoff][e])

			table = odfpy.datatable.DataTable(val)
			table.datasourcehaslabels = "both"
			t = table()
			t.setAttribute('name', "d"+cutoff)
			textdoc.spreadsheet.addElement(t)

	# Table specifique pour un ancetre
	for esp in lstEspeces:
		continue
		val = [["cutoff", "AncGenes", "Blocks", "Genes in blocks", "%Cov", "Min", "25%", "50%", "75%", "N75", "N50", "N25", "Max", "Mean"]]
		for cutoff in allCutoff:
			val.append( [cutoff] + alldata[cutoff][esp][2:-1] )
		
		table = odfpy.datatable.DataTable(val)
		table.datasourcehaslabels = "both"
		t = table()
		t.setAttribute('name', esp)
		textdoc.spreadsheet.addElement(t)


	# Resume final
	val = [["cutoff", "Mean length gain", "N50 gain", "%Cov gain", "Useful Gene Loss"]]
	for cutoff in allCutoff[1:]:
		val.append( [cutoff] + [utils.myMaths.myStats.mean([alldiff[cutoff][e][i] for e in lstEspeces]) for i in [15, 12, 5, 6]] )
	table = odfpy.datatable.DataTable(val)
	table.datasourcehaslabels = "both"
	t = table()
	t.setAttribute('name', "cutoff")
	textdoc.spreadsheet.addElement(t)

	textdoc.save(arguments["outputODS"])

