#!/usr/bin/python
import sys
sys.setrecursionlimit(4000)
import time
import pystp
import os

def readfile(filename):
	f = open(filename, "r")
	tmp = f.readline()
	cur_cost = float(tmp.strip())
	tmp = f.readline()
	N = int(tmp.strip())
	cost_matrix = []
	for i in range(N):
		qq = []
		tmp = f.readline()
		tmp = tmp.strip().split()
		for j in range(N):
			qq.append(float(tmp[j]))
		cost_matrix.append(qq)
	fixed = [{}, {}]
	for _f in range(2):
		tmp = f.readline()
		F = int(tmp.strip())
		for i in range(F):
			tmp = f.readline()
			tmp = tmp.strip().split()
			a = int(tmp[0])
			b = int(tmp[1])
			a,b = sorted([a,b])
			fixed[_f][(b,a)] = float(tmp[2])
	return cur_cost, N, cost_matrix, fixed

def save_solution(filename, solution):
	f = open(filename, "w")
	if solution is None:
		print >>f, "was_optimal"
	else:
		print >>f, "new_solution"
		for x in range(2):
			for i in solution[x]:
				print >>f, i
	f.close()

def preprocess_fixed(fixed):
	sum0 = sum(fixed[0].values())
	sum1 = sum(fixed[1].values())
	delta = sum1 - sum0
	for x in fixed[0].keys():
		fixed[0][x] = 0.0
	for x in fixed[1].keys():
		fixed[1][x] = 0.0
	fixed[1][fixed[1].keys()[0]] = delta
	return fixed, sum0
	

def get_subtours(solution):
	for i in range(N):
		t = 0
		for j in range(i):
			t += solution[i][j]
		for j in range(i+1, N):
			t += solution[j][i]
		assert t == 2

	def get_next(i, prev):
		for j in range(i):
			if solution[i][j] and j != prev:
				return j
		for j in range(i+1, N):
			if solution[j][i] and j != prev:
				return j
		assert False
	used = []
	tours = []
	for i in range(N):
		if i not in used:
			subtour = [i]
			previous = -1
			used.append(i)
			j = get_next(i, previous)
			while j not in used:
				subtour.append(j)
				used.append(j)
				previous = i
				i = j
				j = get_next(i, previous)
			tours.append(subtour)
	tmp = [ x for tour in tours for x in tour ]
	assert len(tmp) == len(solution)
	return tours

def vertex_min_edges(i, cost_matrix, fixed, cnt):
	tmp = []
	tmp_fixed = []
	for j in range(N):
		f = 0
		b,a = sorted([i,j])
		if (a,b) in fixed[0]:
			tmp_fixed.append((fixed[0][(a,b)], (i,j)))
			f += 1
		if (a,b) in fixed[1]:
			tmp_fixed.append((fixed[1][(a,b)], (i,j)))
			f += 1
		if len(fixed) == 3 and (a,b) in fixed[2]:
			tmp_fixed.append((fixed[2][(a,b)], (i,j)))
			f += 1
		assert f <= 2
		if f < 2 and i != j:
			tmp.append((cost_matrix[i][j], (i,j)))
	tmp.sort()
	tmp = tmp_fixed + tmp
	return tmp[:cnt] # fixed + best edges



def get_min_edges(cost_matrix, fixed, cnt):
	edges = []
	for i in range(N):
		tmp = vertex_min_edges(i, cost_matrix, fixed, cnt)
		edges.append(tmp) # save them
	
	cost = sum([cost for l in edges for cost,desc in l])
	return (cost, edges)

def find_irrelevant(cost_matrix, fixed, upper_bound):
	lower_bound, edges = get_min_edges(cost_matrix, fixed, 4)
	# UB cost = max(sum path)
	# LB cost <= sum(path * 2)
	# so adjust UB to match our definition
	upper_bound = upper_bound * 4
	print "min4lb", lower_bound
	print "upperb", upper_bound
	print "delta:", upper_bound - lower_bound
	irrelevant = {}
	for i in range(N):
		for j in range(i):
			if (i,j) in fixed[0] and (i,j) in fixed[1]:
				continue # fixed edge, cannot play with it anyway
			old = edges[i] + edges[j]
			tmp = [ fixed[0],
					fixed[1],
					{(i,j):cost_matrix[i][j]},
				]
			new = vertex_min_edges(i, cost_matrix, tmp, 4)
			new += vertex_min_edges(j, cost_matrix, tmp, 4)
			oldc = sum([cost for cost,desc in old])
			newc = sum([cost for cost,desc in new])
			#print lower_bound, newc, oldc, lower_bound + newc - oldc, upper_bound
			if (lower_bound + newc - oldc - 0.001 > upper_bound):
				irrelevant[(i,j)] = True
				irrelevant[(j,i)] = True
	return irrelevant

def solve2factor_using_subtour_elim(N, fixed, factor4):
	# now we have 4factor, let us use STP
	stp = pystp.Stp()
	two_bit = stp.createType(pystp.TYPE_BITVECTOR, 2)
	
	# constant value
	def Const(value, bits=2):
		return stp.bvConstExprFromLL(bits, value)

	# unsigned less
	def Less(a, b):
		return stp.bvLeExpr(a, b)
	
	# equality
	def Equal(a, b):
		return stp.eqExpr(a, b)

	selected = [None, None]
	print "creating variables"
	for x in range(2):
		# premenna ci je hrana na orientovanej ceste
		selected[x] = []
		for i in range(N):
			selected[x].append([ stp.varExpr("selected_%d_%d_%d" % (x,i,j), two_bit) for j in range(i) ])

	print "assert 0-1 range"
	for x in range(2):
		for i in range(N):
			for j in range(i):
				expr = Less(selected[x][i][j], Const(1))
				stp.assertFormula(expr)

	print "assert 1 for all fixed edges and remove them from 4-factor"
	for x in range(2):
		for i in range(N):
			for j in range(i):
				if (i,j) in fixed[x]:
					expr = Equal(selected[x][i][j], Const(1))
					stp.assertFormula(expr)
					factor4[i][j] = factor4[i][j] - 1 # remove fixed edge from 4-factor
	
	print "assert 0 for all edges not in 4-factor"
	for x in range(2):
		for i in range(N):
			for j in range(i):
				if (i,j) not in fixed[x] and factor4[i][j] == 0:
					formula = Equal(selected[x][i][j], Const(0))
					stp.assertFormula(formula)

	print "assert 2-factors"
	for x in range(2):
		for i in range(N):
			expr = Const(0, 32) # used 32bit variables here
			for j in range(i):
				expr = stp.bvPlusExpr(32, expr, stp.bvConcatExpr(Const(0,30), selected[x][i][j]))
			for j in range(i+1, N):
				expr = stp.bvPlusExpr(32, expr, stp.bvConcatExpr(Const(0, 30), selected[x][j][i]))
			stp.assertFormula(Equal(expr, Const(2, 32)))

	print "assert exclusivity"
	for i in range(N):
		for j in range(i):
			if (i,j) in fixed[0] or (i,j) in fixed[1]:
				continue
			tmp = stp.bvAndExpr(selected[0][i][j], selected[1][i][j])
			expr = Equal(tmp, Const(1))
			stp.assertFormula(stp.notExpr(expr))
	
	print "basic constrains ok, solving"
	
	subtours = []
	subtour_add = []
	while True:
		print "eliminating", len(subtours), "subtours"
		for x in range(2):
			for subtour in subtour_add:
				assert len(subtour) > 2
				# primary condition
				tmp = []
				for i in subtour:
					for j in range(i):
						if j not in subtour:
							tmp.append(selected[x][i][j])
					for j in range(i+1, N):
						if j not in subtour:
							tmp.append(selected[x][j][i])
				assert tmp # wtf tmp is empty?
				expr = Const(0)
				for t in tmp:
					expr = stp.bvOrExpr(expr, t)
				expr = Equal(expr, Const(1)) # or is nonzero
				stp.assertFormula(expr)
		print "querying stp"
		stp.push()
		res = stp.query(stp.falseExpr())
		print res
                if res == True:
                  return None, None
		print "done"
		solution0 = [ [ stp.getCounterExample(selected[0][i][j]).getBVUnsignedLongLong() for j in range(i) ] for i in range(N) ]
		solution1 = [ [ stp.getCounterExample(selected[1][i][j]).getBVUnsignedLongLong() for j in range(i) ] for i in range(N) ]
		stp.pop()
                return solution0, solution1

def load_factor4(factor4file):
  f = open(factor4file, "r")
  out = []
  for l in f:
    if (len(l.strip()) == 0):
      out.append([])
    else:
      out.append([int(x) for x in l.strip().split(' ')]) 
  f.close()
  return out

infile = sys.argv[1]
outfile = sys.argv[2]
factor4file = sys.argv[3]
tempfile = outfile + ".tmp" 
upper_bound, N, cost_matrix,fixed = readfile(infile)
#print N, cost_matrix, fixed
fixed, delta = preprocess_fixed(fixed)
upper_bound -= delta
#irrelevant = find_irrelevant(cost_matrix, fixed, upper_bound)



factor4 = load_factor4(factor4file)
factor41, factor42 = solve2factor_using_subtour_elim(N, fixed, factor4) 
if not factor41:
  save_solution(outfile, None)
  sys.exit(0)

ones = 0
f = open(tempfile, "w")
for i in xrange(0, len(factor41)):
  for j in xrange(0, len(factor41[i])):
    for k in xrange(0, factor41[i][j]):
      print >>f, "%d %d" % (i, j)
    ones += factor41[i][j]
for i in xrange(0, len(factor42)):
  for j in xrange(0, len(factor42[i])):
    for k in xrange(0, factor42[i][j]):
      print >>f, "%d %d" % (i, j)
    ones += factor42[i][j]

print ones
print N

f.close()

retval = os.system("./4factor-to-2paths.bin %s %s %s" % (infile, tempfile, outfile))
print retval


#save_solution(outfile, solution)
