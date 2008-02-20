#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree", "tree.5"], \
	[("IN.CDS",str,""), ("OUT.FASTA-CDS",str,""), ("OUT.FASTA-Prot",str,""), ("OUT.tree",str,"")], \
	"Ecrit les familles" \
)

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree"], buildLinks = False)

CDS = {}
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Chargement des CDS de", esp, "...",
	f = utils.myTools.myOpenFile(options["IN.CDS"] % phylTree.fileName[esp], "r")
	for ligne in f:
		t = ligne.replace('\n','').split('\t')
		CDS[t[0]] = t[1][:3*(len(t[1])/3)] # On tronque aux codons entiers
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
     'GGG': 'G' } #, 'TAA':  '', 'TAG':  '', 'TGA':  '' } # Les codons stop seront des X

f = utils.myTools.myOpenFile(noms_fichiers["tree.5"], "r")
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
		cdsTab[j] = c
		protTab[j] = ''.join([geneticCode.get(c[3*i:3*i+3],'X') for i in xrange(len(c)/3)])
		j += 1
	
	# L'arbre
	tree = utils.myTools.myOpenFile(options["OUT.tree"] % n, "w")
	print >> tree, ligneArbre
	tree.close()
	# Le multi-FASTA des CDS
	fasta = utils.myTools.myOpenFile(options["OUT.FASTA-CDS"] % n, "w")
	print >> fasta, '\n'.join(cdsTab)
	fasta.close()
	# Le multi-FASTA des Proteines
	fasta = utils.myTools.myOpenFile(options["OUT.FASTA-Prot"] % n, "w")
	print >> fasta, '\n'.join(protTab)
	fasta.close()

	print >> sys.stderr, "OK"
f.close()
