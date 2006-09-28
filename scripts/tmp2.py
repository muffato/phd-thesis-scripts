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

def combineSynt(synt1, synt2, seuil):

	assoc = myTools.myCombinator([])
	
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
				assoc.addLink([i,str(j)])
	
	res1 = []
	res2 = []
	for grp in assoc.getGrp():
		res1.append(set(myMaths.flatten([synt1[x] for x in grp if type(x) == int])))
		res2.append(set(myMaths.flatten([synt2[int(x)] for x in grp if type(x) == str])))
	
	print >> sys.stderr, ">", len(res1), sum([len(x) for x in res1])
	return res1

	if seuil == 1:
		return res1

	dic = {}
	for k in range(len(res2)):
		for x in res2[k]:
			dic[x] = k
	for g in res1:
		tmp = set([])
		for x in g:
			if x in dic:
				tmp.update(res2[dic[x]])
		g.update(tmp)
	print >> sys.stderr, ">", len(res1), sum([len(x) for x in res1])

	return res1


# D'abord les syntenies par couples (seuil=1)
newEns = []
esp = geneBank.dicEspeces.values()
for i in range(len(esp)):
	for j in range(len(esp)):
		if j <= i:
			continue
		ens1 = extractSyntenies(esp[i], esp[j], genesAnc)
		ens2 = extractSyntenies(esp[j], esp[i], genesAnc)
		#for a in ens:
		#	if len(a) < options["seuilLongueurMin"]:
		#		continue
		#	print " ".join(myMaths.flatten([genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][u].names for u in a]))

		res = combineSynt(ens1, ens2, 1)
		newEns.append( (esp[i].nom+esp[j].nom,res) )
		#print >> sys.stderr, ">", len(res), sum([len(a) for a in res])
		a = [len(x) for x in res]
		a.sort()
		print >> sys.stderr, a


sys.exit(0)

# Ensuite, on remonte l'arbre
#res = combineSynt(newEns[1][1], newEns[4][1], options["seuilIdentiteMin2"])
#a = [len(x) for x in res]
#print >> sys.stderr, ">>", len(res), sum(a)
#a.sort()
#print >> sys.stderr, a

#newEns2 = []
#for i in range(len(newEns)):
#	for j in range(len(newEns)):
#		if j <= i:
#			continue
#		a = newEns[i]
#		b = newEns[j]
#		res = combineSynt(a[1], b[1], options["seuilIdentiteMin2"])
#		newEns2.append(res)
#		print >> sys.stderr, "%s/%s" % (a[0],b[0]), ">>", len(res), sum([len(x) for x in res])

#sys.exit(0)

#print >> sys.stderr, len(lst), "presque-chromosomes chez l'ancetre"


#print "RESULTAT"
c = 0
n = 0
for aa in res:
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

sys.exit(0)

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
