#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( ["ExpectedAligment", "ComputedAligment"], [], "Compare la sequence reconstruite a la sequence attendue")

expected  = utils.myGenomes.loadFastaFile(noms_fichiers["ExpectedAligment"])
expected["Node1"] = expected["Node1 The root"]
del expected["Node1 The root"]

computed = {}
proba = utils.myTools.defaultdict(dict)
f = utils.myTools.myOpenFile(noms_fichiers["ComputedAligment"], "r")
for l in f:
	t = l.replace('\n','').split('\t')
	if l.startswith("NAME"):
		if t[1] == "NAME_0":
			name = "Node1"
		else:
			name = "Node" + t[1]
	elif l.startswith("SEQ"):
		computed[name] = t[1]
	elif l.startswith("PROBA"):
		base = l[6]
		num = [float(x) for x in t[1].split()]
		proba[base][name] = num
		#print >> sys.stderr, base, name, len(num), num[:50]
f.close()
#print >> sys.stderr, proba.keys()

for anc in computed:
	if len(expected[anc]) != len(computed[anc]):
		#print >> sys.stderr, anc, "DIFFERENT LENGTHS !"
		continue
	nbOK = 0
	pos = 0
	for (i,e) in enumerate(expected[anc]):
		c = computed[anc][i]
		if e == c:
			nbOK += 1
		if c == '-':
			pc = 1
			if e == '-':
				pe = 1
			else:
				pe = 0
		elif c == 'N':
			pc = max([proba[x][anc][pos] for x in proba])
			if e == '-':
				pe = 0
			else:
				pe = proba[e][anc][pos]
			pos += 1
		else:
			pc = proba[c][anc][pos]
			if e == '-':
				pe = 0
			else:
				pe = proba[e][anc][pos]
			pos += 1
		print e, c, pe, pc
	print >> sys.stderr, anc, nbOK, len(expected[anc]),  float(nbOK) / len(expected[anc])
				
		


