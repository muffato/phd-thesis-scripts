#! /usr/bin/python2.4

#
# Charge tous les fichiers d'orthologues lus sur l'entree standard
# Chaque fichier contient les noms des genes sur les colonnes numCol1 et numCol2
# Renvoie la liste des groupes de genes orthologues entre eux
#

import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myTools

(noms_fichiers, options) = myTools.checkArgs([], [("showStats", bool, False)], "Lit sur l'entree standard des familles de genes et les regroupe")

# MAIN #

lst = []
dic = {}

# On scanne toutes les lignes
for l in sys.stdin:
	cc = l.split()
	if len(cc) == 0:
		continue
	c1 = cc[0]
	for c2 in cc[1:]:
		if c1 in dic:
			i = dic[c1]
			if c2 in dic:
			
				# Les deux genes sont deja dans la liste, on les relie
				j = dic[c2]
				
				if i == j:
					continue

				lst[i].extend(lst[j])
				for g in lst[j]:
					dic[g] = i
				lst[j] = []
				
			else:
				# Un seul des deux, on rajoute a cette liste
				lst[i].append( c2 )
				dic[c2] = i
		
		elif c2 in dic:
			# Un seul des deux, on rajoute a cette liste
			j = dic[c2]
			lst[j].append( c1 )
			dic[c1] = j
		
		else:
			# Aucun des deux, on cree un nouveau groupe
			dic[c1] = len(lst)
			dic[c2] = len(lst)
			lst.append( [c1,c2] )

del dic

# On affiche le resultat
for x in lst:
	if len(x) == 0:
		continue
	for g in x:
		print g,
	print
	if options["showStats"]:
		print >> sys.stderr, len(x)

