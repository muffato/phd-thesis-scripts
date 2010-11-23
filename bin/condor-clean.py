#!/usr/bin/env python

# Remove duplicate lines from a condor submission file

import sys

lastValue = {}

for line in sys.stdin:
	if "=" not in line:
		# Don't touch lines without a '='
		print line,
	else:
		(name,_,val) = line.partition("=")
		name = name.strip()
		val = val.strip()
		if (name in lastValue) and (lastValue[name] == val):
			# If the value is the same, do not print the line
			pass
		else:
			# Else, updates the dict
			print line,
			lastValue[name] = val

