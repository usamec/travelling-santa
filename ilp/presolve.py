#!/usr/bin/python
import sys
sys.path.append("/home/ppershing/python/numberjack/local_lib");
sys.setrecursionlimit(4000)
import Numberjack
import time
import pystp

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
					{(i,j):cost_matrix[i][j], (j,i):cost_matrix[i][j]},
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

def get_4factor(N, cost_matrix, fixed, irrelevant, upper_bound, cost_delta):
	print "get_4factor"
	model = Numberjack.Model()
	fix = dict(fixed[0].items() + fixed[1].items())
	vertex_fix = []
	free = {}
	for i in range(N):
		for j in range(i):
			if ((i,j) not in fixed[0] or (i,j) not in fixed[1]):
				if not (i,j) in irrelevant and not (j,i) in irrelevant:
					free[(i,j)] = True
	for i in range(N):
		t = 4;
		for j in range(i):
			if (i,j) in fixed[0]:
				t -= 1
			if (i,j) in fixed[1]:
				t -= 1
		for j in range(i+1, N):
			if (j,i) in fixed[0]:
				t -= 1
			if (j,i) in fixed[1]:
				t -= 1
		vertex_fix.append(t)

	# premenna ci je hrana na orientovanej ceste
	selected = []
	for i in range(N):
		selected.append([ Numberjack.Variable(0, 1) if (i,j) in free else None for j in range(i) ])

	# fix degree
	for i in range(0, N):
		tmp = []
		for j in range(i):
			if (i,j) in free:
				tmp.append(selected[i][j])
		for j in range(i+1, N):
			if (j,i) in free:
				tmp.append(selected[j][i])
		if tmp:
			model.add(Numberjack.Sum(tmp) == vertex_fix[i])
	
	# cost
	a = []
	b = []
	for i in range(N):
		for j in range(i):
			if (i,j) in free:
				a.append(selected[i][j])
				b.append(cost_matrix[i][j])
	cost = Numberjack.Variable(-1e20, 1e20)
	model.add(cost == Numberjack.Sum(a,b))

	objective = Numberjack.Variable(-1e20, 1e20)
	model.add(objective == Numberjack.Sum([cost],[0.5]))
	model.add(Numberjack.Minimize(objective))

	solver = model.load('SCIP')
	solver.setVerbosity(1)
	res = solver.solve()
	assert res == True
	
	lower_bound = objective.get_value()
	lower_bound += sum(fixed[0].values() + fixed[1].values()) / 2
	print ">>>> lower bound", lower_bound + cost_delta
	print ">>>> upper bound", upper_bound + cost_delta
	print ">>>> gap: %.2f " % (upper_bound - lower_bound)
	if lower_bound + 0.01 >= upper_bound:
		print ">>>> was already optimal"
		assert lower_bound - 0.01 < upper_bound
		return None
	solver.delete()
	return True # new solution


infile = sys.argv[1]
outfile = sys.argv[2]
upper_bound, N, cost_matrix,fixed = readfile(infile)
#print N, cost_matrix, fixed
fixed, delta = preprocess_fixed(fixed)
upper_bound -= delta
irrelevant = find_irrelevant(cost_matrix, fixed, upper_bound)



print "irrelevant edges", len(irrelevant)
#print sorted(irrelevant.keys())
factor4 = get_4factor(N, cost_matrix, fixed, irrelevant, upper_bound, delta)
save_solution(outfile, factor4)
