#! /users/ldog/muffato/python -OO

# Fonctions d'appel du voyageur de commerce
# Formats d'entree des donnees:
#  - Fonction de score
# Format de sortie des donnees:
#  - Liste des chemins les plus court


import os
import sys
import random
import utils.myTools


def doConcorde(self, n, func, repeat=1, verbose=False):

	# Pour gerer un resultat de concorde et coorienter plusieurs solutions
	class ConcordeFile:

		# Chargement
		def __init__(self, nom):
			tmp = []
			f = utils.myTools.myOpenFile(nom, 'r')
			for ligne in f:
				for x in ligne.split():
					tmp.append(int(x)-1)
			
			f.close()
			i = tmp.index(-1)
			self.res = tmp[i+1:] + tmp[1:i]
			self.buildIndex()

		# Pour trouver la position d'un element
		def buildIndex(self):
			self.dic = {}
			for (i,x) in enumerate(self.res):
				self.dic[x] = i
		
		# Indice de coorientation de deux solutions
		def isMemeSens(self, other):
			s = 0
			for i in xrange(len(self.res)-1):
				if other.dic[self.res[i]] < other.dic[self.res[i+1]]:
					s += 1
				else:
					s -= 1
			return s>0

		# Inverse le sens de la solution
		def reverse(self):
			self.res.reverse()
			self.buildIndex()
		

	# La matrice des distances intergenes
	print >> sys.stderr, "Ecriture de la matrice ... ",
	filename = "mat%08d" % ((os.getpid() ^ os.getppid() ^ random.randint(1,16777215)) & 16777215)
	f = open(filename, 'w', 65536)

	print >> f, "NAME: CHRANC"
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	print >> f, "0 " * n
	for i in xrange(n):
		for j in xrange(i+1,n):
			print >> f, func(i, j),
		print >> f
	print >> f, "EOF"		
	f.close()
	print >> sys.stderr, "OK"
	
	print >> sys.stderr, "Lancement de concorde ",
	dest = {True:'&2', False:'/dev/null'}[verbose]
	lstTot = []
	for i in xrange(repeat):
		comm = '/users/ldog/muffato/work/scripts/utils/concorde/concorde -x %s >%s' % (filename,dest)
		os.system(comm)
		if utils.myTools.fileAccess(filename + ".sol"):
			lstTot.append(ConcordeFile(filename + ".sol"))
		os.system('rm -f 0%s* %s*' % (filename,filename) )
		sys.stderr.write(".")
	os.system('rm -f *%s*' % filename )

	# On remet chaque liste dans le meme sens que la premiere
	res = []
	for t in lstTot:
		if not t.isMemeSens(lstTot[0]):
			t.reverse()
		res.append(t.res)
	return res


