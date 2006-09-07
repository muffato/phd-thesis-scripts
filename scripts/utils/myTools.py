#! /usr/bin/python2.4


import sys
import os
from bz2 import BZ2File


##################################################################
# Cette classe ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier           #
##################################################################
def myOpenFile(nom):
	nom = nom.replace("~", os.environ['HOME'])
	if nom.endswith("bz2"):
		f = BZ2File(nom, 'r')
	else:
		f = open(nom, 'r')
	return f


#
# Verifie que les arguments obligatoires de la ligne de commande sont presents
# Recupere les valeurs des parametres optionnels
# -options- est un tableau de triplets (nom, constructeur, val_defaut)
#
def checkArgs(args, options, info):

	#
	# Affiche le message d'erreur de mauvais arguments
	#
	def error_usage():
		s = "- ERREUR - Usage : " + sys.argv[0]
		for t in args:
			s += " " + t
		print >> sys.stderr, s
		for t in options:
			if t[1] == bool:
				invite = "+/-"
			else:
				invite = "-"
			print >> sys.stderr, "  ", invite + "%s [%s] (%s)" % t
		if info != "":
			print >> sys.stderr, "\n", info
		sys.exit(1)

	types = dict([ (x[0], x[1]) for x in options ])
	valOpt = dict([ (x[0], x[2]) for x in options ])
	valArg = []
	
	# On scanne les argumetns pour les compter et recuperer les valeurs
	for t in sys.argv[1:]:

		# Un petit peu d'aide
		if t == '-h' or t == '--help':
			error_usage()
			
		# Un argument optionnel
		if t[0] in ['-', '+']:
		
			# Un parametre non bool
			if '=' in t:
				s = t[1:t.index('=')]
				v = t[t.index('=')+1:]
				if not s in valOpt:
					error_usage()
				if types[s] == bool:
					error_usage()
				valOpt[s] = types[s](v)
				
			else:
				s = t[1:]
				if not s in valOpt:
					error_usage()
				if types[s] != bool:
					error_usage()
				# Ici, on affecte False
				valOpt[s] = (t[0] == '+')
		elif os.access(t, os.R_OK):
			valArg.append(t)
		else:
			print >> sys.stderr, "No access to", t
			error_usage()

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) != len(args):
		error_usage()
	
	return (valArg, valOpt)


