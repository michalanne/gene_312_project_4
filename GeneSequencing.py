#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random


# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		if banded:
			score, al1, al2 = self.solveBanded(seq1, seq2, align_length)
		else:
			score, al1, al2 = self.solveUnbanded(seq1, seq2, align_length)
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100
		alignment1 = al2.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = al1.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


	def solveBanded(self, seq1, seq2, align_length):
		k = 7
		# print("seq1: ", seq1, "seq2: ", seq2)
		ysize = len(seq1)
		xsize = len(seq2)
		len1 = len(seq1)
		len2 = len(seq2)
		if align_length < ysize:
			ysize = align_length
		if k < xsize:
			xsize = k # just in case we get a case that is smaller than 7 chars long lol
		prev = []
		array = []
		for ar in range(ysize):
			c = []
			d = []
			for col in range(xsize):
				c.append(0)
				d.append('x')
			array.append(c)
			prev.append(d)
		# E[i,j] = min {diff(i,j)+E[i-1,j-1], 5+E[i,j-1], 5+E[i-1,j]}
		for r in range(ysize):
			array[r][0] = r * 5
		for q in range(xsize):
			array[0][q] = q * 5
		n = 0
		t = 0
		# begin nested for loop that will loop through every single item in array and populate it
		for i in seq1:
			r = 7
			if n >= len1-3:
				r = 4+len1-n
			if n <= 3:
				r = n + 3
			# runs 7 times unless it's the first or last few lines
			for q in range(r):
				d = 1
				if n > 3:
					curr = (n - 3) + t
				else:
					curr = t
				if curr < 0: # error check, but shouldn't need it.
					curr = 0
				if curr >= len(seq2):
					curr = len(seq2) - 1
				if i == seq2[curr]:
					d = -3 # if the items are the same, d = -3, if they're different, d = 1
				if n <= 3:
					diagonal = d+array[n-1][t-1]
					above = 5+array[n-1][t]
					if t == 0:
						side = math.inf # if side does not exist, then it is infinity
					else:
						side = 5+array[n][t-1]
				if n > 3:
					diagonal = d + array[n - 1][t]
					if t == 6:
						above = math.inf # if above node doesn't exist, infinity
					else:
						above = 5 + array[n - 1][t+1]
					if t == 0:
						side = math.inf # if the side does not exist, then it is infinity
					else:
						side = 5 + array[n][t-1]
				if diagonal <= side and diagonal <= above: # compare all of them, set it equal to the smallest
					array[n][q] = diagonal
					prev[n][q] = 's'
				elif side <= above:
					array[n][q] = side
					prev[n][q] = 'i'
				else:
					array[n][q] = above
					prev[n][q] = 'd'
				t += 1
				if t > xsize:
					break
			n += 1
			t = 0
			if n == ysize: # just in case, so it doesn't get too large
				break

		n = ysize - 1
		t = xsize - 1
		al1 = seq2[0:xsize]
		al2 = seq1[0:ysize]
		placeInString1 = xsize
		placeInString2 = ysize
		for cur in range(xsize+ysize):
			if n == 0 and t == 0:
				al1 += seq2[0]
				al2 += seq1[0]
				break
			if prev[n][t] == 's':
				# substitution or equal
				# doesn't mtater if they're equal or substitution, it prints the same thing
				n -= 1
				t -= 1
			elif prev[n][t] == 'i':
				#insertion
				# final_string = my_string[:index] + 'you ' + my_string[index:]
				al1 = al1[:placeInString1] + '-' + al1[placeInString1:]
				t -= 1
			elif prev[n][t] == 'd':
				# deletion
				al2 = al2[:placeInString2] + '-' + al2[placeInString2:]
				n -= 1
			else:
				j = 7
			placeInString1 -= 1
			placeInString2 -= 1

		# a = al1[::-1]
		# a2 = al2[::-1]
		# if seq1 == "exponential" and seq2 == "exponential":
		# 	j = 5/0
		a = al1
		a2 = al2
		a = a[0:100]
		a2 = a2[0:100]

		return array[ysize-1][3], a, a2 # returning finished

	def jB(self, i, j):
		if i <= 3:
			return j
		return j-i+3


	def solveUnbanded(self, seq1, seq2, align_length):
		ysize = len(seq1)
		xsize = len(seq2)
		if align_length < ysize:
			ysize = align_length # make sure you don't assign align_length if it's longer than the length of the word
		if align_length < xsize:
			xsize = align_length
		prev = []
		array = []
		for ar in range(ysize+1): # populate arrays
			c = []
			d = []
			for col in range(xsize+1):
				c.append(0)
				d.append('x')
			array.append(c)
			prev.append(d)
		# E[i,j] = min {diff(i,j)+E[i-1,j-1], 5+E[i,j-1], 5+E[i-1,j]}
		n = 1
		t = 1
		for r in range(ysize+1):
			array[r][0] = r * 5 # populate the first row and the first column so we have something to work with
		for q in range(xsize+1):
			array[0][q] = q * 5
		for i in seq1:
			for j in seq2:
				d = 1
				if i == j:
					d = -3
				# print("d: ", d, " array[n][t-1]: ", array[n][t-1])
				diagonal = d+array[n-1][t-1]
				left = 5+array[n][t-1]
				up = 5+array[n-1][t]
				if left <= diagonal and left <= up:
					array[n][t] = left
					prev[n][t] = 'd'
				elif up <= diagonal:
					array[n][t] = up
					prev[n][t] = 'i'
				else:
					array[n][t] = diagonal
					prev[n][t] = 's'
				# if d+array[n-1][t-1] <= 5+array[n][t-1] and d+array[n-1][t-1] <= 5+array[n-1][t]: # compare above, side, and diagonal
				# 	array[n][t] = d + array[n-1][t-1]
				# 	prev[n][t] = 's'
				# elif 5+array[n][t-1] <= 5+array[n-1][t]:
				# 	array[n][t] = 5+array[n][t-1]
				# 	prev[n][t] = 'i'
				# else:
				# 	array[n][t] = 5+array[n-1][t]
				# 	prev[n][t] = 'd'
				t += 1
				if t > xsize:
					break
			n += 1
			t = 1
			if n > ysize: # just in case
				break
		n = ysize - 1
		t = xsize - 1
		al1 = seq2[0:xsize]
		al2 = seq1[0:ysize]
		placeInString1 = xsize
		placeInString2 = ysize
		# beginning of alignment for loop
		for cur in range(xsize+ysize):
			if prev[n][t] == 's':
				# substitution or equal
				n -= 1
				t -= 1
			elif prev[n][t] == 'i':
				#insertion
				# final_string = my_string[:index] + 'you ' + my_string[index:]
				al1 = al1[:placeInString1] + '-' + al1[placeInString1:]
				n -= 1
			elif prev[n][t] == 'd':
				# deletion
				al2 = al2[:placeInString2] + '-' + al2[placeInString2:]
				t -= 1
			else:
				j = 7
			placeInString1 -= 1
			placeInString2 -= 1
			if n == 0 and t == 0:
				# al1 += seq2[0]
				# al2 += seq1[0]
				break

		# a = al1[::-1]
		# a2 = al2[::-1]
		# if seq1 == "exponential" and seq2 == "exponential":
		# 	j = 5/0
		a = al1
		a2 = al2
		a = a[0:100] # only the first 100 characters
		a2 = a2[0:100]

		return array[ysize][xsize], a, a2
