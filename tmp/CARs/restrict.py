
import sys


data1 = set()

f1 = open(sys.argv[1], "r")
for l in f1:
	data1.update(l[:-1].split()[4:])
f1.close()

f2 = open(sys.argv[2], "r")
for l in f2:
	t = set(l[:-1].split()[4:])
	if len(t.intersection(data1)) > 0:
		print l,
f2.close()

