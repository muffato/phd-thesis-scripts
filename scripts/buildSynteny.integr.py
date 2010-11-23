#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	A partir de diagonales pair-wise, construit des versions integrees qui correspondent a des segments de chromosomes ancestraux
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import utils.myDiags


# Integration de toutes les diagonales d'un ancetre
def getLongestDiags(anc, nbGenesAnc, pairwiseDiags):

	# On construit la table "gene" -> "liste des diagonales qui le contiennent"
	# On s'en sert pour avoir les listes de diagonales chevauchantes
	dic = collections.defaultdict(list)
	combin = utils.myTools.myCombinator()

	for (i,d) in enumerate(pairwiseDiags):
		# Combinaison selon les genes ancestraux
		for (j,_) in d:
			dic[j].append(i)
		combin.addLink([i])

	for j in xrange(nbGenesAnc):
		if j in dic:
			combin.addLink(dic[j])
		else:
			# Envoi des genes singletons
			yield ([j], [1], [])

	del combin.dic
	
	print "new anc", anc
	print "nb genes", nbGenesAnc
	print "in synt", len(dic)

	# Permet d'extraire les composantes connexes
	for g in combin:

		# On casse les carrefours de diagonales les uns apres les autres
		ss = set()
		gr = utils.myDiags.WeightedDiagGraph()
		for i in g:
			ss.update(d for (d,_) in pairwiseDiags[i])
			gr.addDiag(pairwiseDiags[i])
		gr.printIniGraph()
		gr.cleanGraphTopDown(arguments["minimalWeight"])

		# Les diagonales resultat
		tt = set()
		for (res,scores) in gr.getBestDiags():
			strand = [x[1] for x in res]
			res = [x[0] for x in res]
			tt.update(res)
			yield (res, strand, scores)

		# Pour verifier que tous les genes sont bien dans des groupes (eventuellement des singletons)
		assert tt == ss, (ss,tt,ss.difference(tt))


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("minimalWeight",int,1), ("minimalLength",int,2), \
	("IN.projDiags",str,"proj/diags.%s.list.bz2"), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


def do(anc):

	genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	
	print >> sys.stderr, "Chargement des diagonales projetees de %s ..." % anc,
	s = []
	pairwiseDiags = []
	f = utils.myFile.myTSV.reader(arguments["IN.projDiags"] % phylTree.fileName[anc])
	for t in f.csvobject:
		assert len(t) == 10
		assert t[0] == anc
		if (t[2] not in phylTree.listSpecies) or (t[5] not in phylTree.listSpecies):
			continue
		d = zip([int(x) for x in t[8].split()], [int(x) for x in t[9].split()])
		if len(d) >= arguments["minimalLength"]:
			pairwiseDiags.append(d)
			s.append(len(d))
	f.file.close()
	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "OK"

	print >> sys.stderr, "Impression des diagonales ancestrales de %s ..." % anc,
	f = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
	s = []
	singles = 0
	for (da,ds,dw) in getLongestDiags(anc, len(genesAnc.lstGenes[None]), pairwiseDiags):

		if len(da) > 1:
			s.append( len(da) )
		else:
			singles += 1

		res = [anc, len(da), utils.myFile.myTSV.printLine(da," "), utils.myFile.myTSV.printLine(ds, " "), utils.myFile.myTSV.printLine(dw," ")]
		f.csvobject.writerow(res)
	f.file.close()

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "+ %d singletons OK" % singles

# Traitement final
for anc in utils.myDiags.getTargets(phylTree, arguments["target"])[1]:
	do(anc)


