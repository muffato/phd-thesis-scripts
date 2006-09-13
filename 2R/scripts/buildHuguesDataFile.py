#! /usr/bin/python2.4


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

#(noms_fichiers, options) = myTools.checkArgs(["GENOME_CRANIATE", "GENOME_AMNIOTES","GENOME_HUMAIN"], [], "")
(noms_fichiers, options) = myTools.checkArgs(["GENOME_AMNIOTES","GENOME_HUMAIN"], [], "")

#genomeCran = myOrthos.AncestralGenome(noms_fichiers[0], True)
genomeAmn = myOrthos.AncestralGenome(noms_fichiers[0], True)
genomeH = myOrthos.EnsemblGenome(noms_fichiers[1])

nb = 0
for s in sys.stdin:
	t = s.split()
	s = set(t[1:])
	Kcraniate = t[0]
	autres = set([])
	
	Gamn = set([])
	for g in s:
		if g in genomeAmn.dicGenes:
			Gamn.add(genomeAmn.dicGenes[g])
		else:
			autres.add(g)

	if len(Gamn) < 2:
		continue
	
	nb += 1
	print "***** FAMILLE %d *****" % nb
	print "++ %s ++" % " ".join(autres)
	for (Kamn,Iamn) in Gamn:
		s = genomeAmn.lstGenes[Kamn][Iamn][-1]
		Ghum = set([])
		for g in s:
			if g in genomeH.dicGenes:
				Ghum.add(genomeH.dicGenes[g])
		print " ".join(s), "        ", Kcraniate, Kamn, "/".join([str(Khum) for (Khum,_) in Ghum])
	print 
	print 
