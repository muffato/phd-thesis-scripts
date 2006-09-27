#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genesList.conf","genomeOutgroup", "orthologuesList"], \
	[("especesUtilisees",str,"HDWM"), ("seuilLongueurMin",float,0.1), ("seuilIdentiteMin",float,33), ("seuilIdentiteMin2",float,33)], \
	"Reconstruit le genome de l'ancetre de 1 et 2 a partir de l'outgroup et des genes de cet ancetre" \
)


# 1. On lit tous les fichiers
geneBank = myOrthos.GeneBank(noms_fichiers[0], options["especesUtilisees"])
genomeOutgroup = myOrthos.loadGenome(noms_fichiers[-2])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[-1], False)


def extractSyntenies(g1, g2, genesAnc):

	# On construit les couleurs
	res = []
	for c in g1.lstChr:

		tmp = []
		for i in range(len(g1.lstGenes[c])):
		
			tg = myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in g1.lstGenes[c][i].names if g in genesAnc.dicGenes]])
			
			for g in tg:
				if g in g2.dicGenes:
					tmp.append( (g2.dicGenes[g][0],genesAnc.dicGenes[g][1]) )
					break

		last = ""
		curr = set([])
		
		for (col,tg) in tmp:
			if col != last:
				if len(curr) != 0:
					res.append(curr)
				last = col
				curr = set([])
			curr.add(tg)
		res.append(curr)
		
	print >> sys.stderr, g1.nom, g2.nom, len(res), sum([len(a) for a in res])
	return res

class combinator:

	def __init__(self):
		self.dic = {}
		self.grp = []
	
	def addObj(self, obj):
		if len(obj) == 0:
			return
		
		# Le 1er objet de la liste
		a = obj[0]
		if a in self.dic:
			i = self.dic[a]
		else:
			i = len(self.grp)
			self.grp.append([i])
			self.dic[a] = i

		for x in obj[1:]:
			if x in self.dic:
				j = self.dic[x]
				if i == j:
					continue
				tmp = self.grp[j]
				for b in self.grp[j]:
					self.grp[i].append(b)
					self.dic[b] = i
			else:
				self.grp[i].append(x)
				self.dic[x] = i

	def getGrp(self):
		return [g for g in self.grp if len(g) > 0]


def combineSynt(synt1, synt2, seuil):

	assoc = []
	for a in synt1:
		tmp = set([])
		
		for j in range(len(synt2)):
			b = synt2[j]
			z = a & b

			if seuil >= 1:
				test = (len(z) >= seuil)
			else:
				test = (len(z) >= seuil*len(b)) and (len(z) >= seuil*len(a))
			if test:
				tmp.add(j)
		assoc.append(tmp)
	
	ens = range(len(synt1))
	res = []
	
	while len(ens) > 0:
		i = ens.pop()
		a = synt1[i]
		used = True
		while used:
			used = False
			j = 0
			while j < len(ens):
				if len(assoc[i] & assoc[ens[j]]) != 0:
					a.update(synt1[ens[j]])
					del ens[j]
					used = True
				else:
					j += 1
		res.append(a)
	return res

def combineSynt2(synt1, synt2, seuil):

	assoc = combinator()
	
	for i in range(len(synt1)):
		a = synt1[i]
		
		for j in range(len(synt2)):
			b = synt2[j]
			z = a & b

			if seuil >= 1:
				test = (len(z) >= seuil)
			else:
				test = (len(z) >= seuil*len(b)) and (len(z) >= seuil*len(a))
			if test:
				assoc.addObj([i,str(j)])
	
	res = []
	print >> sys.stderr, len(assoc.getGrp())
	for grp in assoc.getGrp():
		for u in [x for x in grp if type(x) == int]:
			print u
		lst = set(myMaths.flatten([synt1[x] for x in grp if type(x) == int]))
		res.append(lst)
	
	return res




newEns = []
for i in range(len(geneBank.dicEspeces.values())):
	for j in range(len(geneBank.dicEspeces.values())):
		if j <= i:
			continue
		ens1 = extractSyntenies(geneBank.dicEspeces.values()[i], geneBank.dicEspeces.values()[j], genesAnc)
		ens2 = extractSyntenies(geneBank.dicEspeces.values()[j], geneBank.dicEspeces.values()[i], genesAnc)
		#for a in ens:
		#	if len(a) < options["seuilLongueurMin"]:
		#		continue
		#	print " ".join(myMaths.flatten([genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][u].names for u in a]))

		res = combineSynt2(ens1, ens2, 1)
		newEns.extend(res)
		print >> sys.stderr, ">", len(res), sum([len(a) for a in res])
		res = combineSynt2(ens2, ens1, 1)
		print >> sys.stderr, ">", len(res), sum([len(a) for a in res])

sys.exit(0)

print >> sys.stderr, len(newEns), sum([len(a) for a in newEns])

lst = []
while len(newEns) > 0:
	a = newEns.pop()
	used = True
	while used:
		used = False
		j = 0
		while j < len(newEns):
			b = newEns[j]
			z = a & b

			#print (100*len(z))/len(a)
			#print (100*len(z))/len(b)
			#print len(z)
			#j += 1
			#continue
			
			if options["seuilIdentiteMin2"] >= 1:
				test = (len(z) >= options["seuilIdentiteMin2"])
			else:
				test = (len(z) >= options["seuilIdentiteMin2"]*len(b))
				test = test or (len(z) >= options["seuilIdentiteMin2"]*len(a))
			if test:
				del newEns[j]
				a.update(b)
				used = True
			else:
				j += 1
	lst.append(a)

#sys.exit(0)

print >> sys.stderr, len(lst), "presque-chromosomes chez l'ancetre"


#print "RESULTAT"
c = 0
n = 0
for aa in lst:
	a = set(aa)
	if len(a) < 10:
		#print >> sys.stderr, len(a),
		continue
	#print len(a)
	#print a
	for i in a:
		#print i
		print chr(97+c), " ".join(genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][i].names)
	n += len(a)
	c += 1

print >> sys.stderr
print >> sys.stderr, n, "genes repartis sur", c, "chromosomes"

sys.exit(1)
res = []
for (l,l1,l2,g1,g2) in lst:

	a = g1 | g2
	if len(a) < 150:
		continue
	res.append( (len(a), l, l1, l2, a) )

print >> sys.stderr, len(res), "presque-chromosomes chez l'ancetre"
res.sort()
res.reverse()
c = 0
nbTMP = 0
for (nb,l,l1,l2,a) in res:
	print >> sys.stderr, l, nb, l1,l2
	for i in a:
		#continue
		print chr(97+c), " ".join(genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][i].names)
	nbTMP += nb
	c += 1
print >> sys.stderr, "Total nb genes:", nbTMP
