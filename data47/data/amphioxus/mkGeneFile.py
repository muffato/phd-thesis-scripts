#! /users/ldog/muffato/python -OO

import sys
import utils.myTools

# scaffold_1  JGI   exon  474   504   .     -     .     name "fgenesh2_pg.scaffold_1000001"; transcriptId 63195

dic = utils.myTools.defaultdict(dict)
f = utils.myTools.myOpenFile(sys.argv[1], 'r')
for l in f:
	t = l.replace('\n', '').split('\t')
	c = t[0]
	x1 = int(t[3])
	x2 = int(t[4])
	strand = t[6]
	name = t[8].split('"')[1]
	if name in dic:
		dic[name][1] = min(dic[name][1], x1)
		dic[name][2] = max(dic[name][2], x1)
	else:
		dic[name] = [c,x1,x2,strand]
f.close()

for (name,(c,x1,x2,strand)) in dic.iteritems():
	print "%s\t%d\t%d\t%d\t%s" % (c,x1,x2,int(strand+"1"),name)

