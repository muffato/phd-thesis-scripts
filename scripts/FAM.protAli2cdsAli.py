#! /users/ldog/muffato/python -OO

import sys
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [], [("range",str,""), ("aligment-FASTA",str,""), ("CDS-FASTA",str,""), ("aligment-output",str,""), ("positions",str,"123")], "Converts a Muscle alignment of protein sequences to the corresponding CDS sequences" )


positions = [int(p)-1 for p in arguments["positions"]]

for treeID in utils.myTools.getRange(arguments["range"]):
	
	print >> sys.stderr, treeID, "...",

	# Chargement des CDS
	CDS = utils.myGenomes.loadFastaFile(arguments["CDS-FASTA"] % treeID)
	# Gestion des CDS non entiers
	for (name,seq) in CDS.iteritems():
		CDS[name] = seq + "NN"

	# Conversion des alignements
	fi = utils.myTools.myOpenFile(arguments["aligment-FASTA"] % treeID, "r")
	fo = utils.myTools.myOpenFile(arguments["aligment-output"] % treeID, "w")
	seqFiltered = ""
	for ligne in fi:
		ligne = ligne.replace('\n', '')
		if ligne[0] == ">":
			print >> fo, seqFiltered
			print >> fo, ligne
			name = ligne[1:]
			seqFiltered = ""
			pos = 0
		else:
			# Parcours des acides amines
			for aa in ligne:
				# Recuperation du codon associe
				if aa == '-':
					codon = '---'
				else:
					codon = CDS[name][pos:pos+3]
					pos += 3
				# On imprime les positions voulues
				seqFiltered += ''.join([codon[p] for p in positions])
	print >> fo, seqFiltered
	fi.close()
	fo.close()
	print >> sys.stderr, "OK"

