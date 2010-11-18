#!/usr/bin/env python2

import sys

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [], [("range",str,""), ("aligment-FASTA",str,""), ("CDS-FASTA",str,""), ("aligment-output",str,""), ("positions",str,"123")], "Converts a Muscle alignment of protein sequences to the corresponding CDS sequences" )


positions = [int(p)-1 for p in arguments["positions"]]

for treeID in utils.myTools.getRange(arguments["range"]):
	
	print >> sys.stderr, treeID, "...",

	# Chargement des CDS
	CDS = utils.myGenomes.myFASTA.loadFile(arguments["CDS-FASTA"] % treeID)
	# Gestion des CDS non entiers
	for (name,seq) in CDS.iteritems():
		CDS[name] = seq + "NN"

	# Conversion des alignements
	fi = utils.myFile.openFile(arguments["aligment-FASTA"] % treeID, "r")
	fo = utils.myFile.openFile(arguments["aligment-output"] % treeID, "w")
	for ligne in fi:
		ligne = ligne.replace('\n', '')
		if ligne[0] == ">":
			print >> fo, ligne
			name = ligne[1:]
			pos = 0
		else:
			seqFiltered = []
			# Parcours des acides amines
			for aa in ligne:
				# Recuperation du codon associe
				if aa == '-':
					codon = '---'
				else:
					codon = CDS[name][pos:pos+3]
					pos += 3
				# On imprime les positions voulues
				seqFiltered.extend( codon[p] for p in positions )
			print >> fo, ''.join(seqFiltered)
	fi.close()
	fo.close()
	print >> sys.stderr, "OK"

