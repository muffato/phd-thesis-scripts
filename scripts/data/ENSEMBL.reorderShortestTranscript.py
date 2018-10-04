#!/usr/bin/env python2

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("genesFile",file), ("transcriptsCoords",file)], \
	[("useShortestTranscript",bool,True), ("sortOn5",bool,True), ("printTranscriptCoords",bool,True)], \
	"Cree une liste ordonnee des genes en tenant compte du plus petit transcrit" \
)

genome = utils.myGenomes.Genome(arguments["genesFile"])

# Chargement de la liste des transcrits
lstTrans = collections.defaultdict(list)
f = utils.myFile.myTSV.reader(arguments["transcriptsCoords"])
for l in f.csvobject:
	lstTrans[l[0]].append( (int(l[2]),int(l[3]),l[1]) )
f.file.close()

for chrom in genome.lstGenes:

	# Creation de la liste a trier
	tmp = []
	for gene in genome.lstGenes[chrom]:
		if arguments["useShortestTranscript"]:
			assert gene.names[0] in lstTrans, gene.names[0]
			l = [(y-x,x,y,t) for (x,y,t) in lstTrans[gene.names[0]]]
			best = min(l)
			tmp.append(best[1:3] + (gene,best[3]))
		else:
			tmp.append((gene.beginning,gene.end,gene,None))

	# Tri selon le 5' ou le 3'
	if arguments["sortOn5"]:
		tmp.sort()
	else:
		import operator
		tmp.sort(key=operator.itemgetter(1))

	# Affichage
	for (x,y,gene,name) in tmp:
		res = [gene.chromosome, x, y, gene.strand, " ".join(gene.names)]
		if arguments["printTranscriptCoords"]:
			res[1] = x
			res[2] = y
		if name is not None:
			res.append(name)
		print utils.myFile.myTSV.printLine(res)

