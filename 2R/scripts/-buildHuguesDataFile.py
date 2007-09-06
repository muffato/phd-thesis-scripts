#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

(noms_fichiers, options) = myTools.checkArgs(["genomeAmniotes","genomeHumain"], [], \
	"Affiche toutes les familles de genes craniates separement avec les chromsoomes amniotes et humains sur lesquels ces genes se trouvent")

genomeAmn = myOrthos.AncestralGenome(noms_fichiers[0], True, False)
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
		s = genomeAmn.lstGenes[Kamn][Iamn].names
		Ghum = set([])
		for g in s:
			if g in genomeH.dicGenes:
				Ghum.add(genomeH.dicGenes[g])
		print " ".join(s), "        ", Kcraniate, Kamn, "/".join([str(Khum) for (Khum,_) in Ghum])
	print 
	print 
