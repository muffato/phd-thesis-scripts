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
	["geneList.conf","genomeOutgroup", "orthologuesList"], \
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
	
	res = res1 #[x for x in res1 if len(x) >= 10]
	print >> sys.stderr, ">", len(res), sum([len(x) for x in res])
	return res

	#if seuil == 1:
	#	return res1

	#dic = {}
	#for k in range(len(res2)):
	#	for x in res2[k]:
	#		dic[x] = k
	#for g in res1:
	#	tmp = set([])
	#	for x in g:
	#		if x in dic:
	#			tmp.update(res2[dic[x]])
	#	g.update(tmp)
	#print >> sys.stderr, ">", len(res1), sum([len(x) for x in res1])

	#return res1

def mixGenomes(g1, g2):
	ens1 = extractSyntenies(g1, g2, genesAnc)
	ens2 = extractSyntenies(g2, g1, genesAnc)
	return combineSynt(ens1, ens2, 1)


def buildAncestrGenome(genomesList, lstGenes, assoc):

	print >> sys.stderr, "1.",
	tabE = dict([(e,{}) for e in genomesList])
	for i in xrange(len(lstGenes)):
		g = lstGenes[i]
		for s in g.names:
			for e in genomesList:
				if s in genomesList[e].dicGenes:
					(c,_) = genomesList[e].dicGenes[s]
					tabE[e][i] = c
					break

	print >> sys.stderr, "2.",
	for (i,j) in myTools.myMatrixIterator(len(lstGenes), len(lstGenes), myTools.myMatrixIterator.StrictUpperMatrix):
		# Si les deux genes sont deja rassembles,
		#   on peut passer a la suite
		if assoc.dic[i] == assoc.dic[j]:
			continue
		nbOK = set([])
		nbNO = set([])
		for e in genomesList:
			if (i not in tabE[e]) or (j not in tabE[e]):
				continue
			if tabE[e][i] == tabE[e][j]:
				nbOK.add(e)
			else:
				nbNO.add(e)
		
		# Les trois branches: DW/HM/OC
		if ('H' in nbOK or 'M' in nbOK) and ('D' in nbOK or 'W' in nbOK) and ('O' in nbOK or 'C' in nbOK):
			assoc.addLink([i,j])
		elif 'H' in nbOK and 'M' in nbOK and 'D' in nbOK and 'W' in nbOK:
			assoc.addLink([i,j])

	print >> sys.stderr, "OK"
	#return assoc



def translateGenome(genome, genesAnc, default):

	# On construit les couleurs
	newGen = {}
	for c in genome.lstChr:

		newChr = []
		for gene in genome.lstGenes[c]:
			
			for g in gene.names:
				if g in genesAnc.dicGenes:
					newChr.append( genesAnc.dicGenes[g][1] )
					break
			else:
				newChr.append(default)
		newGen[c] = newChr
	return newGen


genomeH = translateGenome(geneBank.dicEspeces['H'], genesAnc, 'H')
genomeM = translateGenome(geneBank.dicEspeces['M'], genesAnc, 'M')
genomeD = translateGenome(geneBank.dicEspeces['D'], genesAnc, 'D')
genomeW = translateGenome(geneBank.dicEspeces['W'], genesAnc, 'W')
genomeO = translateGenome(geneBank.dicEspeces['O'], genesAnc, 'O')
genomeC = translateGenome(geneBank.dicEspeces['C'], genesAnc, 'C')


def useDiags(genome1, genome2, assoc):

	for c1 in genome1:
		for c2 in genome2:
			diags = myMaths.extractDiags(genome1[c1], genome2[c2], False, 0)
			for d in diags:
				if len(d) >= 2:
					assoc.addLink(d)
	print >> sys.stderr, "*", len(assoc.getGrp())
	
#assoc = myTools.myCombinator([])
lstGenes = genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr]
assoc = myTools.myCombinator([set([u]) for u in xrange(len(lstGenes))])
useDiags(genomeH, genomeD, assoc)
useDiags(genomeH, genomeW, assoc)
useDiags(genomeM, genomeD, assoc)
useDiags(genomeM, genomeW, assoc)
useDiags(genomeH, genomeO, assoc)
useDiags(genomeM, genomeO, assoc)
useDiags(genomeD, genomeO, assoc)
useDiags(genomeW, genomeO, assoc)
useDiags(genomeH, genomeC, assoc)
useDiags(genomeM, genomeC, assoc)
useDiags(genomeD, genomeC, assoc)
useDiags(genomeW, genomeC, assoc)
buildAncestrGenome(geneBank.dicEspeces, lstGenes, assoc)

#assoc = buildAncestrGenome(geneBank.dicEspeces.values(), genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr])
#res = buildAncestrGenome(geneBank.dicEspeces.values(), genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr])


res = assoc.getGrp()

#resDW = mixGenomes(geneBank.dicEspeces['D'], geneBank.dicEspeces['W'])
#resHM = mixGenomes(geneBank.dicEspeces['H'], geneBank.dicEspeces['M'])
#resCO = mixGenomes(geneBank.dicEspeces['C'], geneBank.dicEspeces['O'])
#resCT = mixGenomes(geneBank.dicEspeces['C'], genomeOutgroup)
#resHC = mixGenomes(geneBank.dicEspeces['H'], geneBank.dicEspeces['C'])

#resMam = combineSynt(resHM, resDW, options["seuilIdentiteMin2"])
#resAmn = combineSynt(resMam, resCO, options["seuilIdentiteMin2"])
#resAmn2 = combineSynt(resHM, resCT, options["seuilIdentiteMin2"])

#res = resHC

# D'abord les syntenies par couples (seuil=1)
#newEns = []
esp = geneBank.dicEspeces.values() # + [genomeOutgroup]
for i in range(len(esp)):
	continue
	for j in range(len(esp)):
		if j <= i:
			continue
		res = mixGenomes(esp[i], esp[j])
		#for a in ens:
		#	if len(a) < options["seuilLongueurMin"]:
		#		continue
		#	print " ".join(myMaths.flatten([genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][u].names for u in a]))
		#newEns.append( (esp[i].nom+esp[j].nom,res) )
		for x in res:
			assoc.addLink([u for u in x])
		#newEns.append(res)
		#a = [len(x) for x in res]
		#a.sort()
		#print >> sys.stderr, a


#res = combineSynt(newEns[0], newEns[1], 1)

#res = assoc.getGrp()

#sys.exit(0)

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

#res = newEns[0][1]
c = 0
n = 0
for a in res:
	if len(a) < 2:
		#print >> sys.stderr, len(a),
		continue
	n += len(a)
	c += 1
	print len(a)
	#print a
	continue
	for i in a:
		#print i
		print chr(97+c), " ".join(genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][i].names)

print >> sys.stderr
print >> sys.stderr, n, "genes repartis sur", c, "chromosomes"

