#! /users/ldog/muffato/python -OO

import os
import sys
import math
import time
#import numpy
import random
import zipfile
#import operator
import utils.myGenomes
#import utils.myTools
import utils.myMaths
#import utils.myDiags
#import utils.myPsOutput
import utils.myPhylTree
#import utils.walktrap
#from collections import defaultdict

for i in xrange(10):
	print os.isatty(sys.stdin.fileno())
	time.sleep(1)

sys.exit(0)
print utils.myTools.myIterator.buildSubsets(range(4), 1)
print utils.myTools.myIterator.buildSubsets(range(4), 2)
print utils.myTools.myIterator.buildSubsets(range(4), 3)
print utils.myTools.myIterator.buildSubsets(range(4), 4)

print utils.myTools.myIterator.buildSubsets(range(4), -1)
#print utils.myTools.myIterator.buildSubsets(range(5), 4)

sys.exit(0)
phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], False)
phylTree2 = utils.myPhylTree2.PhylogeneticTree(sys.argv[1])

random.seed(0)

#values = dict.fromkeys(phylTree.listSpecies, 0)

s1 = []
t1 = 0
s2 = []
t2 = 0
print >> sys.stderr, len(phylTree.listSpecies)
for _ in xrange(10000):
	values = {}
	for esp in random.sample(phylTree.listSpecies, 7):
		values[esp] = 1
	t1 -= time.clock()
	s1.append(phylTree.calcWeightedValue(values, -5, None, None))
	t1 += time.clock()
	t2 -= time.clock()
	s2.append(phylTree2.calcWeightedValue(values, -5, None))
	t2 += time.clock()

print >> sys.stderr, all([all(s1[i]==s2[i]) for i in xrange(len(s1))]), t1, t2

sys.exit(0)


#for i in xrange(100):
#phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], True)
phylTree2 = utils.myPhylTree2.PhylogeneticTree(sys.argv[1])

#print phylTree.dicParents == phylTree2.dicParents
#print phylTree.dicLinks == phylTree2.dicLinks

#print phylTree.dicParents["Theria"]
#print phylTree.dicLinks["Theria"]


sys.exit(0)
print phylTree.listAncestr

for e in phylTree.listSpecies:
	values = dict.fromkeys(phylTree.listSpecies, 1)
	values[e] = 1/math.e
	r = phylTree.calcWeightedValue(values, -5, None, None)
	#print e, " ".join([ "%.2g" % (r[phylTree.indNames[a]]) for a in phylTree.listAncestr ])
	print e, "%g" % (1/r[phylTree.indNames["Boreoeutheria"]])

sys.exit(0)



f = utils.myGenomes.loadFastaFile("/workspace2/muffato/ancSequences47/tmp/cds.fa")

for (name,seq) in f.iteritems():
	print name
	print len(seq)
	print seq

sys.exit(0)
f = utils.myTools.myOpenFile("~/family.1.zip/ancSequence.txt", "r")
nbL = 0
for l in f:
	print l,
	nbL += 1

print >> sys.stderr, nbL

sys.exit(0)

# ESSAIS POUR randomSlice



length = 1000.
nb = [0] * int(1.5*length)

for _ in xrange(1000000):
	#r = math.pow(2., random.uniform(-1., 1.))
	r = math.pow(2., random.vonmisesvariate(0, 1) / math.pi)
	nb[int(length*(r-.5))] += 1

for (i,v) in enumerate(nb):
	print (i/length)+.5, v

sys.exit(0)

values = {}
#values["Tetraodon nigroviridis"] = 1
values["Gallus gallus"] = 1
values["Homo sapiens"] = 1

print >> sys.stderr, phylTree.allNames
print >> sys.stderr, phylTree.calcWeightedValue(values, -5, "Amniota", None)

sys.exit(0)


phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], True)
values = {}
for espR in phylTree.listSpecies:
	for esp in phylTree.listSpecies:
		values[esp] = 0
	values[espR] = 1
	for anc in phylTree.listAncestr:
		r = phylTree.calcWeightedValue(values, -5, None, anc)
		if r[1] != anc:
			print >> sys.stderr, "!", r, espR, anc
		r = r[2]
		comm = phylTree.dicParents[espR][anc]
		d = 2*phylTree.ages[comm] - phylTree.ages[anc]
		if comm == anc:
			tmp = espR
			d = 0
			while tmp != anc:
				(par,dist) = phylTree.parent[tmp]
				if len(phylTree.items[par]) > 1:
					dist *= len(phylTree.items[par])
				d += dist
				tmp = par
			print espR, anc, r, 1./d, r*d
				

sys.exit(0)


phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], False)
print phylTree.items
sys.exit(0)

# Write

#f = zipfile.ZipFile("tata.zip", "r", zipfile.ZIP_STORED)
f = zipfile.ZipFile("toto.zip", "w", zipfile.ZIP_STORED)

s = "bonjour, jeune camarade !\nComment vas-tu ?\n"

zinfo = zipfile.ZipInfo(filename="f1l3n4m3", date_time=time.localtime(time.time()))
zinfo.compress_type =  zipfile.ZIP_DEFLATED
zinfo.external_attr = 2175008768
f.writestr( zinfo, s*100000)

f.writestr( "f1l3n4m3", s*100)

f.close()
sys.exit(0)

inf = f.getinfo("scripts/FamFetch/Help")
#print inf
#print dir(inf)
for a in dir(inf):
	print a, getattr(inf,a)
#print f.namelist()
#print f.infolist()
#print f.read("scripts")

inf = f.getinfo("scripts/buildPreDupGenome.py")
for a in dir(inf):
	print a, getattr(inf,a)

f.close()
sys.exit(0)

# ESSAIS POUR randomSlice

length = 1000
mean = .2
nb = [0] * (length+600)
x0 = 2*mean - 1
if x0 > 0:
	x0 = -x0

for _ in xrange(1000000):
	# r entre 0 et 1
	r = random.vonmisesvariate(0,2) / (2*math.pi) + .5
	r = abs(x0 + r*(1-x0))
	#r = x0 + r*(1-x0)
	if mean > 0.5:
		r = 1 - r
	nb[int(length*r+600)] += 1

for (i,v) in enumerate(nb):
	print i-600, v

sys.exit(0)

length = 10000

# 5000 - 7500 - 2500 / 27"
def v1():
	x1 = random.randint(0, length-1)
	x2 = random.randint(x1, length-1)
	return (x1,x2)


# 36"
# 1 = 3400 -> 6600 = 3200 
# 2 = 3900 -> 6100 = 2200
# 3 = 4200 -> 5800 = 1600
def v2():
	l = int(abs(random.vonmisesvariate(0, 4)) * length / math.pi)
	x1 = random.randint(0, length-1-l)
	return (x1,x1+l)


# 2500 - 7500 - 5000 / 27"
def v3():
	l = random.randint(0, length-1)
	x1 = random.randint(0, length-1-l)
	return (x1,x1+l)


# 3333 - 6666 - 3333 / 30"
def v4():
	x = random.sample(xrange(length), 2)
	return (min(x),max(x))


# 3333 - 6666 - 3333 / 30"
def v5():
	x1 = random.randint(0, length-1)
	x2 = random.randint(0, length-1)
	return (min(x1,x2),max(x1,x2))

res1 = []
res2 = []
resL = []
for i in xrange(1000000):
	(x1,x2) = v2()
	res1.append(x1)
	res2.append(x2)
	resL.append(x2-x1+1)
print utils.myMaths.myStats(res1)
print utils.myMaths.myStats(res2)
print utils.myMaths.myStats(resL)


sys.exit(0)
for _ in xrange(00):
	phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], False)
#dsi = dict.__setitem__

for _ in xrange(000):
	r = phylTree.newCommonNamesMapperInstance()
	for (a,t) in phylTree.officialName.iteritems():
		#try:
			dict.__setitem__(r, t, a)
			#dsi(r, t, a)
			#r.setdefault(t, a)
			#r[t] = a
			#phylTree.species[t]
		#except KeyError:
		#	pass

phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1], False)
values = {}
#values["ENSGACG00000000274"] = 0
values["Tetraodon nigroviridis"] = 0
values["Danio rerio"] = 1

print phylTree.allNames
print phylTree.calcWeightedValue(values, -5, "Clupeocephala", "Stickleback")
sys.exit(0)

print
print phylTree.allNames
print
print phylTree.indNames
print
print phylTree.items
print
print phylTree.root
print
print phylTree.parent
print
print

print values

print
print
print phylTree.calcWeightedValue(values, -5, "Clupeocephala", None)
print phylTree.calcWeightedValue(values, -5, "Clupeocephala", "Clupeocephala")

print phylTree.allNames
print
print phylTree.indNames
print
print phylTree.items
print
print phylTree.root
print
print phylTree.parent
print
print

print values

print phylTree.calcWeightedValue(values, -5, "Clupeocephala", "Clupeocephala")
print phylTree.calcWeightedValue(values, -5, "Clupeocephala", None)

sys.exit(0)


n = 4000
#l = [[0] * n for _ in xrange(n)]
#l = [0] * n
#l = numpy.array(l)
l = numpy.zeros( (n,n) )
print l

for i in xrange(n):
	#m[i] = i
	#continue
	t = l[i]
	for j in xrange(n):
		#m[i,j] = i+j
		#l[i][j] = i+j
		t[j] = i+j

sys.exit(0)
import numpy
for _ in xrange(10000):
	for _ in xrange(1000):
		import numpy
		#pass
		#import numpy

#phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

#print len(phylTree.allNames)

#print 1
#print phylTree.parent["Human"]

#print 2
#print phylTree.parent.get("Human", "*NON*")

#print 3
#print phylTree.parent["Homo sapiens"]

sys.exit(0)
phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

tree = utils.myTree.NewickTree("(W09G3.3:0.912111,((((SINFRUG00000163035:0.064359,GSTENG00018392001:0.055428):0.033866,ENSORLG00000010372:0.159161):0.906856,ENSGALG00000015890:0.591867):0.383249,(ENSCSAVG00000000038:0.333015,ENSCING00000004007:0.200122):0.386033):0.0) ;")
print tree.data
print tree.items
print tree.root
print tree.species
print


tree = utils.myTree.NewickTree("(W09G3.3 25.71:0.6125, ((ENSCSAVG00000000038 48.37:0.2672, ENSCING00000004007 36.11:4.7298)62.52:4.415183, (((SINFRUG00000163035 70.98:0.0768, GSTENG00018392001 74.30:0.1940)66.21:0.495336, ENSORLG00000010372 76.38:0.3924)93.88:7.135469, ENSGALG00000015890 29.44:0.0372)31.38:0.089361)27.19:2.718108);")
print tree.data
print tree.items
print tree.root
print

values = {}
for x in tree.allNames:
	if ' ' in x:
		gc = float(x[x.index(' ') + 1:])
		values[x] = gc
print values

print tree.calcWeightedValue(values)

sys.exit(0)
#import utils.psyco
#utils.psyco.full()

class RTE:
	pass
	#def __init__(self):
	#	self.allNames = []
	#	self.items = {}
	#	self.root = None

self = RTE()
self.allNames = ["A","B","C","D","1","2","3","4","5"]
self.items["A"] = [("1",1),("2",1)]
self.items["B"] = [("A",1),("3",2)]
self.items["C"] = [("B",1),("D",2)]
self.items["D"] = [("4",1),("5",1)]
self.root = "C"

values = {}
values["1"] = 1
values["2"] = 2
values["3"] = 4
values["4"] = 1
values["5"] = 2

print utils.myTree.calcWeightedValue(self, values)

sys.exit(0)
f1 = utils.myTools.myOpenFile(sys.argv[1], 'r')
f2 = utils.myTools.myOpenFile(sys.argv[2], 'r')
lx1 = []
lx2 = []
for l1 in f1:
	l2 = f2.next()
	try:
		x1 = float(l1)
		x2 = float(l2)
		lx1.append(x1)
		lx2.append(x2)
		print x1, x2
	except ValueError:
		pass

print >> sys.stderr, len(lx1)
print >> sys.stderr, utils.myMaths.correlation(lx1, lx2)

sys.exit(0)

x = range(100)
y = x
#y = range(1000, 0, -1)

print x
print y
print utils.myMaths.correlation(x, y)


sys.exit(0)

f = utils.myTools.myOpenFile(sys.argv[1], 'r')
nb = 1
for i in xrange(42370):
	ligne = f.readline()
	arb = f.readline()
	ng = len(ligne.split())
	if ng >= 3:
		toto = open("/workspace/muffato/EVAL_NDGC/%d/tree.txt" % nb, "w")
		print >> toto, arb,
		print nb
		toto.close()
		nb += 1

sys.exit(0)
genome1 = utils.myTools.myOpenFile(sys.argv[1], 'r')
genome2 = utils.myTools.myOpenFile(sys.argv[2], 'r')

for (i,g) in enumerate(genome1):
	g1 = frozenset(g.split())
	g2 = frozenset(genome2.readline().split())
	if g1 != g2:
		print "CONTENT-MISMATCH", i, g1, g2
	continue
	pos = genome2.getPosition(g.names)
	if len(pos) != 1:
		print "LENGTH-MISMATCH", i, g, pos
	else:
		(_,j) = pos.pop()
		if i != j:
			print "INDEX-MISMATCH", i, g, j



sys.exit(0)

seq = sys.stdin.readline()[:-1]
for l in [1,2,3,4,5,6,7]:
	count = utils.myTools.defaultdict(int)
	for x in xrange(len(seq)-l+1):
		count[seq[x:x+l]] += 1
	print utils.myMaths.myStats(count.values()), 1000000 / (4**l)
sys.exit(0)

bases = "ACGT"
seq = ""
nb = 0
#f = open("/dev/urandom", "rb", 0)
#f = random.Random()
#seq = []
#for b in bases:
#	seq.extend( [b] * 250000)
#random.shuffle(seq)
#seq = "".join(seq)
for _ in xrange(1000000):
	break
	seq += f.choice(bases)
	nb += 1
while nb < 1000000:
	input = sys.stdin.read(500)
	#break
	#input = os.urandom(250000)
	for c in input:
		c = ord(c)
		for _ in xrange(4):
			x = c & 3
			c >>= 2
			seq += bases[x]
		nb += 4
#f.close()
print seq[:1000000]
sys.exit(0)

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


