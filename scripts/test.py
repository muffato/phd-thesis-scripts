#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
#import os
import sys
#import math
#import time
#import numpy
#import random
#import operator
#import utils.myGenomes
import utils.myTools
import utils.myMaths
#import utils.myDiags
#import utils.myPsOutput
#import utils.myPhylTree
#import utils.walktrap
#from collections import defaultdict

#import utils.psyco
#utils.psyco.full()

for l in sys.stdin:
	x = int(l.split()[0])
	#print l[:-1], "%e" % utils.myMaths.proba(1./12., x, 487)
	#print l[:-1], "%e" % sum([utils.myMaths.proba(1./12., y, 487) for y in xrange(x,488) ])
	print l[:-1], "%e" % sum([utils.myMaths.proba(1./12., y, 487) for y in xrange(x+1) ])

sys.exit(0)



keys = range(int(sys.argv[1]))
nb = int(sys.argv[2])


def getSubsets1(keys, nb):
	res = set()
	lst = [keys] * nb
	for t in utils.myTools.myIterator.tupleOnManyLists(*lst):
		tt = frozenset(t)
		if tt in res:
			continue
		res.add(tt)
	return res
def buildSubsets(lst, n):
	l = len(lst)
	mem = {}

	def rec(i, n):
		if (i,n) in mem:
			return mem[(i,n)]
		
		if i >= l-n:
			res = [lst[i:]]
		elif n <= 1:
			res = [[x] for x in lst[i:]]
		else:
			ref = [lst[i]]
			res = rec(i+1, n)
			for x in rec(i+1, n-1):
				res.append(ref + x)

		mem[(i,n)] = res
		return res

	return rec(0, n)
	
#print len(getSubsets1(keys, nb))
print len(utils.myTools.myIterator.buildSubsets(keys, nb))
sys.exit(0)

f = utils.myTools.myOpenFile(sys.argv[1], "r")
nb = {}
for l in f:
	t = l[:-1].split('\t')
	t = "|".join(t[4:])
	r = []
	for x in t.split('|'):
		if x.startswith("Homo sapiens"):
			r.append( x.split('/')[1] )
	print " ".join(sorted(r))
sys.exit(0)
f = utils.myTools.myOpenFile(sys.argv[1], "r")
nb = {}
for l in f:
	t = l[:-1].split('\t')
	t = "|".join(t[4:])
	e = [x.split('/')[0] for x in t.split('|')]
	for x in set(e):
		if x not in nb:
			nb[x] = defaultdict(int)
		nb[x][e.count(x)] += 1
for e in nb:
	print e
	print dict(nb[e])
sys.exit(0)

f = utils.myTools.myOpenFile(sys.argv[1], "r")
bad = set()
for l in f:
	bad.update(l.split())
f.close()
print len(bad)

f = utils.myTools.myOpenFile(sys.argv[2], "r")
for l in f:
	t = l.split()
	x = bad.intersection(t)
	if len(x) > 0:
		print x

sys.exit(0)

p = [1,2,5,3,1]

randPick = utils.myMaths.randomValue(p)
nb = [0] * len(p)
for i in xrange(10000000):
	nb[randPick.getRandomPos()] += 1

print nb
print [nb[i]/p[i] for i in xrange(len(p))]

sys.exit(0)

phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

def fileName(anc):
	if anc in phylTree.listAncestr:
		return "/users/ldog/muffato/work/ancestralGenomes43/3.sorted/Genome.%s.bz2" % phylTree.fileName[anc]
	else:
		return "/users/ldog/muffato/work/data43/genes/genes.%s.list.bz2" % phylTree.fileName[anc]

for anc in phylTree.listAncestr:

	for (esp,_) in phylTree.items[anc]:
		res = [anc, esp]
		(stdin,stdout,stderr) = os.popen3("/users/ldog/muffato/work/scripts/printDiags.py +psyco -sameStrand -fusionThreshold=0 %s %s +includeScaffolds +includeRandoms -orthologuesList=/users/ldog/muffato/work/data43/ancGenes/ancGenes.%s.list.bz2" % (fileName(anc),fileName(esp),phylTree.fileName[anc]))
		stdin.close()
		for l in stdout:
			pass
		res = None
		for l in stderr:
			if l.startswith("Extraction des diagonales"):
				t = l.split()
				res = [t[7], t[5].split('/')[1], t[6].split('/')[1]]
				t = t[8][1:-1].replace('/',' ').replace('-',' ').split()
				res += [t[0],t[2]]
				break
		#if res == None:
		#	continue
		res = [anc, esp, phylTree.ages[anc], phylTree.ages[esp]] + res
		print "\t".join([str(x) for x in res])
		stdout.close()
		stderr.close()


sys.exit(0)


def fileName(anc):
	if anc in phylTree.listAncestr:
		return "/users/ldog/muffato/work/ancestralGenomes43/1.ini/Genome.%s.bz2" % phylTree.fileName[anc]
	else:
		return "/users/ldog/muffato/work/data43/genes/genes.%s.list.bz2" % phylTree.fileName[anc]

for anc in phylTree.listAncestr:

	for (esp,_) in phylTree.items[anc]:
		res = [anc, esp]
		(stdin,stdout,stderr) = os.popen3("/users/ldog/muffato/work/scripts/printOrthologousChr.py +psyco %s %s +includeScaffolds +includeRandoms" % (fileName(anc),fileName(esp)))
		stdin.close()
		nbC = "0"
		nbF = "0"
		for l in stderr:
			if "Nb points de cassures" in l:
				nbC = l.split()[4]
			elif "Nb points de fusions" in l:
				nbF = l.split()[4]
		res = [anc, esp, phylTree.ages[anc], phylTree.ages[esp], nbC, nbF]
		print "\t".join([str(x) for x in res])
		stdout.close()
		stderr.close()


sys.exit(0)


for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.listSpecies):

	(stdin,stdout) = os.popen2("/users/ldog/muffato/work/scripts/searchConservedSynteny.py %s +psyco -species=%s,%s -genesFile=data46/genes/genes.%%s.list.bz2 -ancGenesFile=data46/ancGenes/ancGenes.%%s.list.bz2 | awk -f /users/ldog/muffato/work/scripts/searchConservedSynteny.mkOrientationClass.awk | sort | awk '$1==$2' | uniq -c" % (sys.argv[1],e1.replace(' ','^'),e2.replace(' ','^')))

	stdin.close()
	resC = resD = resU = 0
	for ligne in stdout:
		t = ligne.split()
		if t[1] == 'C':
			resC = float(t[0])
		elif t[1] == 'D':
			resD = float(t[0])
		elif t[1] == 'U':
			resU = float(t[0])
	s = (resC+resD+resU)/100.
	if s == 0:
		continue
	print phylTree.ages[phylTree.dicParents[e1][e2]], phylTree.fileName[e1], phylTree.fileName[e2], phylTree.fileName[phylTree.dicParents[e1][e2]], resC/s, resD/s, resU/s


sys.exit(0)
#print len(utils.myTools.myIterator.buildSubsets(range(int(sys.argv[1])), int(sys.argv[2])))

f = open(sys.argv[1], "r")
for ligne in f.xreadlines():
	print "*%s*" % ligne
	break
f.close()

f = open(sys.argv[1], "r")
for ligne in f.readlines():
	print "*%s*" % ligne
	break
f.close()

f = open(sys.argv[1], "r")
for ligne in f:
	print "*%s*" % ligne
	break
f.close()

sys.exit(0)

comb = utils.myTools.myCombinator()

random.seed(int(sys.argv[1]))

ll = []
for i in xrange(100000):
	
	l = [random.randint(1,10000000) for i in xrange(random.randint(2,100))]

	#ll.extend(l)
	comb.addLink(l)

#print comb.getNbGrp()
#print len(utils.myMaths.flatten(ll))
#print len(ll)
#print utils.myMaths.myStats(ll)

#time.sleep(10)

sys.exit(0)


n = 2000
#x = [ range(5000) for x in xrange(5000)]
x = numpy.zeros( (n,n) )

#for (i,j) in utils.myTools.myIterator.tupleOnStrictUpperList(range(n)):
#	print x[i][j]
for i in xrange(n):
	tmp = x[i]
	for j in xrange(i):
		print tmp[j]
		#print x[(i,j)]


sys.exit(0)

g = utils.myGenomes2.Genome(sys.argv[1])
print g.lstChr
print g.lstScaff
print g.lstRand
print len(list(g))

sys.exit(0)

phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

print phylTree.dicParents['Boreoeutheria']['Oryzias latipes']

sys.exit(0)


s = set()
r = range(1000)
#for t in utils.myTools.myMatrixIterator(r, None, utils.myTools.myMatrixIterator.WholeMatrix):
#for t in utils.myTools.combination(r,r,r+r):
for t in utils.myTools.myIterator.tupleOnWholeList(r):
	s.add(t)

print s.pop()

sys.exit(0)

utils.myTools.checkArgs([],[],"")

lst = [tuple(xrange(100)) for i in xrange(200000)]
#lst = [list(xrange(100)) for i in xrange(200000)]

time.sleep(100)
sys.exit(0)

print lst
print utils.myMaths.myStats(lst)

def randomPlace():
	r = random.randint(0, sum(lst)-1)
	for (c,l) in enumerate(lst):
		if r < l:
			return (c, r)
		r -= l


for i in xrange(1000):
	(x1,_) = randomPlace()
	(x2,_) = randomPlace()
	#(x1,x2) = random.sample(range(20), 2)
	if lst[x1] >= 5:
		tmp = lst[x1] / 5
		lst[x1] -= tmp
		lst[x2] += tmp
		
print lst
print utils.myMaths.myStats(lst)

#time.sleep(100)
sys.exit(0)

phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])


def getSumOfBranchLength(node):
	
	s = 0
	for (f,x) in phylTree.items.get(node,[]):
		s += x
		s += getSumOfBranchLength(f)
	return s
		

for anc in phylTree.listAncestr:

	out = []
	tmp = anc
	while tmp in phylTree.parent:
		par = phylTree.parent[tmp]
		out.extend([e for e in phylTree.branches[par] if e != tmp])
		tmp = par

	#nbEsp = [len(phylTree.species[x]) for (x,age) in phylTree.items[anc]]
	nbEspEq = [float(age+getSumOfBranchLength(x))/phylTree.ages[anc] for (x,age) in phylTree.items[anc]]
	nbOutEq = [float(getSumOfBranchLength(x))/phylTree.ages[x] for x in out if x in phylTree.listAncestr] + [1. for x in out if x in phylTree.listSpecies]
	nbEspEq.append(sum(nbOutEq))
	print anc, nbEspEq
	print anc, sum([x*y for (x,y) in utils.myTools.myMatrixIterator(nbEspEq, None, utils.myTools.myMatrixIterator.StrictUpperMatrix)]),
	print math.sqrt(sum([(x-sum(nbEspEq)/len(nbEspEq))**2 for x in nbEspEq]))/len(nbEspEq)

sys.exit(0)

def linear(a, b, ratio):
	return int(a + ratio*(b-a))

def linear3((ra,ga,ba), (rb,gb,bb), ratio):
	return (linear(ra,rb,ratio), linear(ga,gb,ratio), linear(ba,bb,ratio))

def printLine(start, end, nbsteps, y):
	for i in xrange(nbsteps):
		col = linear3(start, end, i/(nbsteps-1.))
		utils.myPsOutput.drawBox(i*21./nbsteps, y, 21/nbsteps, 1, col, col)

chrC = {}
chrC["1"] = (153,102,0)
chrC["2"] = (102,102,0)
chrC["3"] = (153,153,30)
chrC["4"] = (204,0,0)
chrC["5"] = (255,0,0)
chrC["6"] = (255,0,204)
chrC["7"] = (255,204,204)
chrC["8"] = (255,153,0)
chrC["9"] = (255,204,0)
chrC["10"] = (255,255,0)
chrC["11"] = (204,255,0)
chrC["12"] = (0,255,0)
chrC["13"] = (53,128,0)
chrC["14"] = (0,0,204)
chrC["15"] = (102,153,255)
chrC["16"] = (153,204,255)
chrC["17"] = (0,255,255)
chrC["18"] = (204,255,255)
chrC["19"] = (153,0,204)
chrC["20"] = (204,51,255)
chrC["21"] = (204,153,255)
chrC["22"] = (102,102,102)
chrC["23"] = (153,153,153)
chrC["24"] = (204,204,204)


utils.myPsOutput.printPsHeader()


#printLine( (10,10,100), (230,230,255), 3, 1)
#printLine( (10,10,100), (10,10,255), 3, 2)

#printLine( (10,100,10), (230,255,230), 3, 4)
#printLine( (100,10,10), (255,230,230), 3, 3)

#printLine( (10,10,255), (230,230,255), 3, 5)
#printLine( (10,255,10), (230,255,230), 3, 6)
#printLine( (255,10,10), (255,230,230), 3, 7)
#printLine( (10,255,255), (230,255,255), 3, 8)
##printLine( (255,10,255), (255,230,255), 3, 9)
#printLine( (255,255,10), (255,255,230), 3, 10)

#for i in xrange(1,25+1):
#	utils.myPsOutput.drawBox(i*.5, 15, .5, 1, i, i)

#for i in chrC:
#	utils.myPsOutput.drawBox(int(i)*.5, 16, .5, 1, chrC[i], chrC[i])

c5 = []
c5 += ["red4", "coral4", "firebrick", "red2", "coral2", "darkorange", "gold"]
c5 += ["yellow", "khaki1","wheat1","peachpuff","lightsalmon", "hotpink2", "magenta2", "darkorchid2", "purple2", "darkorchid4"]
c5 += ["blue2", "royalblue2", "blue4"]
c5 += ["turquoise4", "darkseagreen4", "chartreuse4","mediumaquamarine",(121,204,61), "chartreuse2", "olivedrab2", "darkolivegreen1"]
c5 += ["darkseagreen1", "paleturquoise2", "lightblue", "skyblue1", "turquoise2", "lavender","thistle2"]
c5 += [(204,204,153),"lightgoldenrod3","ivory2","honeydew3","slategray"]
for i in xrange(len(c5)):
	utils.myPsOutput.drawBox(i*.5, 21, .5, 1, c5[i], c5[i])




utils.myPsOutput.printPsFooter()

sys.exit(0)

def func(x):
	return x*x

for i in xrange(100):

	test = [random.random() for j in xrange(100000)]
	#todo = [x*x for x in test]
	#todo = map(lambda x: x*x, test)
	#todo = map(func, test)
	todo = range(len(test))
	for j in xrange(len(test)):
		todo[j] = test[j]*test[j]


sys.exit(0)



phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

for anc in phylTree.listAncestr:
	
	nbO = len(phylTree.outgroupSpecies[anc])
	f = [len(x) for x in phylTree.branchesSpecies[anc]]
	s = nbO * sum(f)
	for (n1,n2) in utils.myTools.myMatrixIterator(f, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		s += n1*n2

	print anc, s


sys.exit(0)
nbb = 10

x = "1"
x = x.zfill(int(sys.argv[1]))
xx = [int(x[i:i+nbb],2) for i in xrange(0,len(x)+1-nbb,nbb)]

y = "1"
y = y.zfill(int(sys.argv[2]))
yy = [int(y[i:i+nbb],2) for i in xrange(0,len(y)+1-nbb,nbb)]

dic = {}
for i in xrange(2**nbb):
	dic[i] = {}
	for j in xrange(2**nbb):
		n = i ^ j
		tmp = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111)
		dic[i][j] = ((tmp + (tmp >> 3)) & 030707070707) % 63
tyy = [dic[i] for i in yy]

for i in xrange(0, len(xx)-len(yy)):
	s = sum([tyy[j][xx[i+j]] for j in xrange(len(yy))])
	print i, 2*s-len(y)

sys.exit(0)

n = len(y)
while len(x) >= n:

	s = len([None for j in xrange(-1,-n-1,-1) if x[j] == y[j]])
	print len(x)-n, 2*s-n
	x.pop()

sys.exit(0)

#genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.Gallus.gallus.list.bz2")
#genome = utils.myGenomes.EnsemblGenome(sys.argv[1])
phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])
phylTree.initCalcDist("Boreoeutheria", True)
#phylTree.buildPhylLinks()
#print phylTree.dicLinks["Homo sapiens"]
#print phylTree.dicLinks["Boreoeutheria"]
#print phylTree.dicParents["Homo sapiens"]
#print phylTree.dicParents["Boreoeutheria"]

print phylTree.tmpItems

sys.exit(0)

for e in phylTree.listSpecies:

	genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.%s.list.bz2" % phylTree.fileName[e])
	
	s = 0
	lst = {}
	for c in genome.lstChr + genome.lstScaff:
		s += len(genome.lstGenes[c])
		lst[len(genome.lstGenes[c])] = lst.get(len(genome.lstGenes[c]), 0) + 1
		

	lst2 = lst.items()
	lst2.sort(reverse=True)

	p = s
	nbc = 0
	for (l,n) in lst2:
		nbc += n
		p -= l
		if p*100 < s*1:
			break

	print nbc, "\t", e

sys.exit(0)

genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.Gallus.gallus.list.bz2")
#for g in genome.lstGenes[4]:
#	print g.beginning, g.strand, g.names

genesAnc = utils.myGenomes.EnsemblOrthosListGenome("~/work/data43/orthologs/orthos.Gallus.gallus.Homo.sapiens.list.bz2", ancFilter = ["Amniota"])
(c,i) = genesAnc.dicGenes["ENSGALG00000008006"]
print genesAnc.lstGenes[c][i].names
#print genesAnc.dicGenes["ENSGALG00000007993"].names
#print genesAnc.dicGenes["ENSGALG00000008006"].names
#print genesAnc.dicGenes["ENSG00000077279"].names
#print genesAnc.dicGenes["ENSG00000077274"].names

sys.exit(0)

utils.myPsOutput.printPsHeader()
utils.myPsOutput.drawLine(10, 10, 1, 1, utils.myPsOutput.getColor(1, "black"))
utils.myPsOutput.drawLine(10, 10, -1, 1, utils.myPsOutput.getColor(2, "black"))
utils.myPsOutput.drawLine(10, 10, -1, -1, utils.myPsOutput.getColor(3, "black"))
utils.myPsOutput.drawLine(10, 10, 1, -1, utils.myPsOutput.getColor(4, "black"))
utils.myPsOutput.printPsFooter()


sys.exit(0)
for i in xrange(8):
	for j in xrange(3):
		s = -1
		if (i & (1 << j)) > 1:
			s = 1
		print i, j, j, s, "gene.%d.%d" % (i,j)



sys.exit(0)
genesAnc = utils.myGenomes.AncestralGenome(sys.argv[1])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)

for i in xrange(1000):
	print "A", " ".join(lstGenesAnc[i].names)
sys.exit(0)

w = utils.myCommunities2.WalktrapWrapper()

#f = open('/workspace/muffato/proteines/reducedGraphs/graph.0', 'r')
f = open('/users/ldog/muffato/work/graph.3109', 'r')

w.updateFromFile(f)
w.doWalktrap()

#print w.res
#sys.exit(0)

for (comp,best,d) in w.res:

	print >> sys.stderr, "comp connexe"
	for (scale,rel) in best:
		print >> sys.stderr, "scale", scale
		print >> sys.stderr, "relevance", rel
		print >> sys.stderr, d.cut(scale)
sys.exit(0)


for i in xrange(100000):
	dic[random.random()] = random.random()
	dic[random.random()] = str(random.random())
	dic[random.random()] = 1000000*random.random()
	dic[str(random.random())] = random.random()
	dic[str(random.random())] = str(random.random())
	dic[str(random.random())] = 1000000*random.random()
	dic[1000000*random.random()] = random.random()
	dic[1000000*random.random()] = str(random.random())
	dic[1000000*random.random()] = 1000000*random.random()

for (k,v) in dic.iteritems():
	print k, v

sys.exit(0)


def countAltern(lst):

	# La liste des chromosomes de l'alternance
	#lst = [x[eD] for (_,_,x) in lstDCS.__reversed__()]

	# Compte le nombre d'occurrences de c dans la liste courante
	def countChr(c):
		nb = 0
		for x in lst:
			if c not in x:
				break
			nb += 1
		return nb

	# Le compte final
	count = defaultdict(int)
	# La derniere ligne lue
	last = defaultdict(int)
	# On parcourt la liste
	while len(lst) > 0:
		curr = lst.pop()
		for x in curr:
			# Les alternances sont mesurees entre deux positions consecutives
			for y in last:
				if y == x:
					continue
				count[(x,y)] += (countChr(x)+1) * last[y]
				count[(y,x)] = count[(x,y)]
			# Et aussi entre les paralogues
			for y in curr:
				if y >= x:
					continue
				count[(x,y)] += 1
				count[(y,x)] = count[(x,y)]

		# On met a jour last
		for y in last:
			if y not in curr:
				last[y] = 0
		for x in curr:
			last[x] += 1

	return count

#print countAltern([[5],[13],[5],[5],[5,13],[5,13],[5],[13]])
print countAltern([[5],[6],[5],[5],[5,13],[7,13],[13],[5],[1],[5],[5],[13],[13],[1,5]])


sys.exit(0)


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf" ], [], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])


print phylTree.parent['Homo sapiens']
print phylTree.parent['Human']
print phylTree.species['Boreoeutheria']
print phylTree.species['Glires']
print phylTree.branchesSpecies['Vertebreta']

print phylTree.fileName


