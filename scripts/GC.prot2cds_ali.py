#! /users/ldog/muffato/python -OO

import sys
import utils.myTools
import utils.myGenomes

(_,options) = utils.myTools.checkArgs( [], [("start",int,0), ("end",int,0), ("aligment-FASTA",str,""), ("CDS-FASTA",str,""), ("aligment-output",str,""), ("positions",str,"123")], "Converts a Muscle alignment of protein sequences to the corresponding CDS sequences" )

for treeID in xrange(options["start"], options["end"]+1):
	
	print >> sys.stderr, treeID, "...",

	# Chargement des CDS
	CDS = utils.myGenomes.loadFastaFile(options["CDS-FASTA"] % treeID)

	# Conversion des alignements
	fi = utils.myTools.myOpenFile(options["aligment-FASTA"] % treeID, "r")
	fo = utils.myTools.myOpenFile(options["aligment-output"] % treeID, "w")
	seqFiltered = ""
	for ligne in fi:
		if ligne[0] == ">":
			print >> fo, seqFiltered
			print >> fo, ligne,
			name = ligne[1:-1]
			seqFiltered = ""
			pos = 0
		else:
			# Parcours des acides amines
			for aa in ligne[:-1]:
				# Recuperation du codon associe
				if aa == '-':
					codon = '---'
				else:
					codon = (CDS[name][pos:pos+3] + "---")[:3] # Astuce pour les CDS avec une taille non multiple de 3
					pos += 3
				# On imprime les positions voulues
				for p in options["positions"]:
					seqFiltered += codon[int(p)-1]
	print >> fo, seqFiltered
	fi.close()
	fo.close()
	print >> sys.stderr, "OK"

