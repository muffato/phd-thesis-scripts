#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import os
import sys
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree",file), ("tree.5",file)], \
	[("IN.CDS",str,""), ("OUT.FASTA-CDS",str,""), ("OUT.FASTA-Prot",str,""), ("OUT.tree",str,"")], \
	"Ecrit les familles" \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])

CDS = {}
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Chargement des CDS de", esp, "...",
	f = utils.myTools.myOpenFile(arguments["IN.CDS"] % phylTree.fileName[esp], "r")
	for ligne in f:
		t = ligne.replace('\n','').split('\t')
		CDS[t[0]] = t[1] + "NN" # Pour tenir compte des CDS non entiers
	f.close()
	print >> sys.stderr, "OK"


geneticCode = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*' }

f = utils.myTools.myOpenFile(arguments["tree.5"], "r")
for (i,ligneEspeces) in enumerate(f):
	n = i+1
	print >> sys.stderr, "Ligne", n,
	ligneArbre = f.next().replace('\n', '')

	# Generation des donnees FASTA prot & cds
	esp = [(g,CDS[g]) for g in ligneEspeces.split()]

	cdsTab = [0] * (2*len(esp))
	protTab = [0] * (2*len(esp))
	j = 0
	for (g,c) in esp:
		cdsTab[j] = protTab[j] = '>' + g
		j += 1
		length = len(c)/3
		cdsTab[j] = c[:3*length]
		protTab[j] = ''.join([geneticCode.get(c[3*i:3*i+3],'X') for i in xrange(length)])
		# Le * final est enleve car non gere par muscle
		if protTab[j][-1] == '*':
			protTab[j] = protTab[j][:-1]
		# Il faut pour la meme raison remplacer les * par des X
		protTab[j] = protTab[j].replace('*', 'X')
		j += 1

	try:
		# L'arbre
		tree = utils.myTools.myOpenFile(arguments["OUT.tree"] % n, "w")
		print >> tree, ligneArbre
		tree.close()
		# Le multi-FASTA des CDS
		fasta = utils.myTools.myOpenFile(arguments["OUT.FASTA-CDS"] % n, "w")
		print >> fasta, '\n'.join(cdsTab)
		fasta.close()
		# Le multi-FASTA des Proteines
		fasta = utils.myTools.myOpenFile(arguments["OUT.FASTA-Prot"] % n, "w")
		print >> fasta, '\n'.join(protTab)
		fasta.close()
	except IOError:
		print >> sys.stderr, "IOError",

	print >> sys.stderr, "OK"
f.close()
