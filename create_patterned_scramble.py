#!/usr/bin/python

import random
import sys

def nature(aa):
	if aa == 'A' : return 'O'
	if aa == 'C' : return 'C'
	if aa == 'D' : return 'B'
	if aa == 'E' : return 'B'
	if aa == 'F' : return 'O'
	if aa == 'G' : return 'G'
	if aa == 'H' : return 'B'
	if aa == 'I' : return 'O'
	if aa == 'K' : return 'B'
	if aa == 'L' : return 'O'
	if aa == 'M' : return 'O'
	if aa == 'N' : return 'B'
	if aa == 'P' : return 'P'
	if aa == 'Q' : return 'B'
	if aa == 'R' : return 'B'
	if aa == 'S' : return 'B'
	if aa == 'T' : return 'B'
	if aa == 'V' : return 'O'
	if aa == 'W' : return 'O'
	if aa == 'Y' : return 'O'

seq = sys.argv[1]

# Store identities:
O = []
B = []
for aa in seq:
	if nature(aa) == 'O':
		O.append(aa)
	elif nature(aa) == 'B':
		B.append(aa)
	else:
		continue
# Create scramble:
scramble = ''
random.shuffle(O)
random.shuffle(B)
for aa in seq:
	if nature(aa) == 'O':
		scramble += O.pop()
	elif nature(aa) == 'B':
		scramble += B.pop()
	else:
		scramble += aa

print scramble
