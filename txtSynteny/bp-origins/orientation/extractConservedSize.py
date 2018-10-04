#!/usr/bin/env python2

import sys
import math
import collections

import utils.myMaths

last = None

def pr():
	l1 = utils.myMaths.myStats.getValue(lst["True"], 50)
	l2 = utils.myMaths.myStats.getValue(lst["False"], 50)
	l1 = utils.myMaths.myStats.mean(lst["True"])
	l2 = utils.myMaths.myStats.mean(lst["False"])
	if (l1 is not None) and (l2 is not None):
		print " ".join(last), l1, l2, math.log(l1/float(l2), 2)

for l in sys.stdin:
	t = l[:-1].split()
	if t[:5] != last:
		if last is not None:
			pr()
		lst = collections.defaultdict(list)
		last = t[:5]
	if int(t[5]) > 0:
		lst[t[6]].append(int(t[5]))
pr()


