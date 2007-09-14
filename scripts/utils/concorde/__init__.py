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


class ConcordeLauncher:

	def __init__(self):

		self.filename = "mat%08d" % ((os.getpid() ^ os.getppid() ^ random.randint(1,16777215)) & 16777215)
	
	def doConcorde(self, n, func, repeat=1, verbose=False):

		# La matrice des distances intergenes
		print >> sys.stderr, "Ecriture de la matrice ... ",
		f = open(self.filename, 'w', 65536)

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
			comm = '/users/ldog/muffato/work/scripts/utils/concorde/concorde -x %s >%s' % (self.filename,dest)
			os.system(comm)
			if os.access(self.filename + ".sol", os.R_OK):
				lstTot.append(ConcordeFile(self.filename + ".sol"))
			os.system('rm -f 0%s* %s*' % (self.filename,self.filename) )
			sys.stderr.write(".")
		os.system('rm -f *%s*' % self.filename )

		# On remet chaque liste dans le meme sens que la premiere
		res = []
		for t in lstTot:
			if not t.isMemeSens(lstTot[0]):
				t.reverse()
			res.append(t.res)
		return res

	
#####################################################
# Cette classe gere un fichier resultat de Concorde #
#####################################################
class ConcordeFile:

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
			

	def buildIndex(self):
		self.dic = {}
		for (i,x) in enumerate(self.res):
			self.dic[x] = i
	
	def isMemeSens(self, other):
		s = 0
		for i in xrange(len(self.res)-1):
			if other.dic[self.res[i]] < other.dic[self.res[i+1]]:
				s += 1
			else:
				s -= 1
		return s>0

	def reverse(self):
		self.res.reverse()
		self.buildIndex()
	


