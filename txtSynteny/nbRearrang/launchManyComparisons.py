#! /users/ldog/muffato/python

__doc__ = """
	A partir d'un genome de reference et de 3 genomes voisins, extrait la liste des rearrangements specifiques de chaque branche
"""

import sys
import itertools
import subprocess
import collections

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("confFile",file), ("modernGenomesFiles",str), ("ancestralGenomesFiles",str)], \
	[("exec",file,"./extractAllAdjacencies.py")], \
	__doc__ \
)

# Ancetre ou espece moderne
def filename(s):
	if s[0] == "*":
		return arguments["ancestralGenomesFiles"] % s[1:]
	else:
		return arguments["modernGenomesFiles"] % s

# Association de deux noms de genes a la liste des intervalles qu'ils englobent
def getInterv(genome, (g1, g2)):
	if g1 not in genome.dicGenes:
		assert "/" in g1
		return []
	if g2 not in genome.dicGenes:
		assert "/" in g2
		return []
	(c1,i1) = genome.dicGenes[g1]
	(c2,i2) = genome.dicGenes[g2]
	assert c1 == c2
	if i1 < i2:
		return [(c1,i,i+1) for i in xrange(i1,i2)]
	else:
		return [(c1,i,i+1) for i in xrange(i2,i1)]


dic12 = {}
dic21 = {}

def setIntoDict(dic, key, val):

	if key in dic:
		assert dic[key] == val, (key, dic[key], val)

def doCmp(speciesList, lossPattern, gainPattern):

	nesp = len(speciesList)
	ngains = len(gains)
	nloss = len(loss)
	iLost = speciesList.index(esp1)
	iGain = speciesList.index(esp2)

	for (d,l) in itertools.product(parDup, parLength):
		args = [arguments["exec"], "+psyco", "-minLength=%d" % l, d + "eliminateDup"] + [filename(x) for x in speciesList]
		p = subprocess.Popen(args, stdout=subprocess.PIPE)
		for x in p.stdout:
			t = x.split()
			if all(x in y for (x,y) in zip(t, gainPattern)):
				gains.append( getInterv(genome2, (t[nesp+1+iGain], t[2*nesp+1+iGain])) )
				setIntoDict(dic21, t[nesp+1+iGain], t[nesp+1+iLost])
				setIntoDict(dic21, t[2*nesp+1+iGain], t[2*nesp+1+iLost])
				#dic21[t[nesp+1+iGain]] = t[nesp+1+iLost]
				#dic21[t[2*nesp+1+iGain]] = t[2*nesp+1+iLost]
			elif all(x in y for (x,y) in zip(t, lossPattern)):
				loss.append( getInterv(genome1, (t[nesp+1+iLost], t[2*nesp+1+iLost])) )
				setIntoDict(dic12, t[nesp+1+iLost], t[nesp+1+iGain])
				setIntoDict(dic12, t[2*nesp+1+iLost], t[2*nesp+1+iGain])
				#dic12[t[nesp+1+iLost]] = t[nesp+1+iGain]
				#dic12[t[2*nesp+1+iLost]] = t[2*nesp+1+iGain]
		p.wait()
		print >> sys.stderr, "gains:", len(gains)-ngains, "loss:", len(loss)-nloss

# Combinaison des intervalles et selection des intervalles les plus souvent vus
def selectBest(lstInt, txt, genome, dic):

	# Nom du gene dans les deux genomes (si disponible)
	def getName(g):
		s = g.names[0]
		if s not in dic:
			return s
		t = dic[s]
		if t == s:
			return s + "|="
		return s + "|" + t

	print >> sys.stderr, txt, len(lstInt)
	count = collections.defaultdict(int)
	comb = utils.myTools.myCombinator()
	for interv in lstInt:
		comb.addLink(interv)
		for x in interv:
			count[x] += 1
	print >> sys.stderr, txt, len(comb.dic), len(comb.grp), len(list(comb))
	for (i,grp) in enumerate(comb):
		grp.sort()
		scores = [(count[x], genome.lstGenes[x[0]][x[1]], genome.lstGenes[x[0]][x[2]]) for x in grp]
		bestscore = max(scores)[0]
		extr = [(g1,g2) for (s,g1,g2) in scores if s == bestscore]
		print utils.myFile.myTSV.printLine([
			txt, "chr:%s" % grp[0][0],
			getName(extr[0][0]), getName(extr[-1][1]), "score=%d" % bestscore, "len=%d" % len(extr), "orient=%d/%d" % (extr[0][0].strand,extr[-1][1].strand),
			"{ " + getName(scores[0][1]) + "".join(" [%d] %s" % (s,getName(g2)) for (s,g1,g2) in scores) + " }" ])

# Lecture du fichier de conf
loss = []
gains = []
f = utils.myFile.openFile(arguments["confFile"], "r")
for l in f:

	if l.startswith(">"):
		if (len(loss) > 0) or (len(gains) > 0):
			# Enregistrement
			selectBest(loss, "LOST", genome1, dic12)
			selectBest(gains, "GAIN", genome2, dic21)

		# Nouvelles especes de references et nouveaux parametres
		(_, esp1, esp2, parDup, parLength) = l[:-1].split("\t")
		genome1 = utils.myGenomes.Genome(filename(esp1))
		genome2 = utils.myGenomes.Genome(filename(esp2))
		parLength = [int(x) for x in parLength.split()]
		loss = []
		gains = []


	elif l.startswith("#"):
		# Commentaire
		pass

	else:
		# Nouvelle liste d'especes a comparer et motifs de reconnaissance des rearrangements
		(speciesList, lossPattern, gainPattern) = [x.split() for x in l[:-1].split("\t")]
		assert len(lossPattern) == len(speciesList)
		assert len(gainPattern) == len(speciesList)
		doCmp(speciesList, lossPattern, gainPattern)

f.close()
if (len(loss) > 0) or (len(gains) > 0):
	# Enregistrement
	selectBest(loss, "LOST", genome1, dic12)
	selectBest(gains, "GAIN", genome2, dic21)

