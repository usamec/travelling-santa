#!/usr/bin/python
import sys
sys.path.append("/home/ppershing/python/numberjack/local_lib");
sys.setrecursionlimit(4000)
import Numberjack
import time

assert len(sys.argv) == 3

fin_filename = sys.argv[1]
fout_filename = sys.argv[2]

fin = open(fin_filename, "r")
tmp = fin.readline()
tmp = fin.readline();
N = int(tmp.strip())
print "Problem dimension: ", N
external_edge_cost = []
for i in range(N):
	tmp = fin.readline()
	tmp = tmp.strip().split(' ')
	assert len(tmp) == 4
	external_edge_cost.append([])
	for j in range(4):
		external_edge_cost[i].append(float(tmp[j]))


normal_edge_cost = []
for i in range(N):
	tmp = fin.readline()
	tmp = tmp.strip().split(' ')
	assert(len(tmp) == N)
	normal_edge_cost.append([])
	for j in range(N):
		normal_edge_cost[i].append(float(tmp[j]))

fin.close()

for i in range(N):
	for j in range(N):
		assert(normal_edge_cost[i][j] == normal_edge_cost[j][i])

def get_min_edges(normal_edge_cost, external_edge_cost, cnt_edges, cnt_skip):
	edges = []
	flat = []
	for i in range(N):
		tmp = []
		for x in external_edge_cost[i]:
			tmp.append((x, (i, 'external')))
		tmp += [(normal_edge_cost[i][j], (i,j)) for j in range(i)]
		tmp += [(normal_edge_cost[i][j], (i,j)) for j in range(i+1,N)]
		tmp.sort()
		tmp = tmp[:cnt_edges] # best edges
		edges.append(tmp) # save them
		flat += tmp 
	flat.sort()
	flat.reverse()
	for i in range(cnt_skip):
		edge = flat[i]
		f = edge[1][0] # from
		edges[f].remove(edge)
	
	cost = sum([cost for l in edges for cost,desc in l])
	return (cost, edges, flat[:cnt_skip])


def solve(normal_edge_cost, external_edge_cost, cnt_skip):
	min10cost, min10edges, _ = get_min_edges(normal_edge_cost, external_edge_cost, 10, cnt_skip)
	forbidden = {}
	for i in range(N):
		if (len(min10edges[i]) != 10):
			continue
		maxc = max(min10edges[i])[0]

		for j in range(i):
			if normal_edge_cost[i][j] > maxc:
				forbidden[(i,j)] = True
				forbidden[(j,i)] = True
	#print "getting UB"
	cost10 = solve_using_lp(normal_edge_cost, external_edge_cost, cnt_skip, forbidden)
	
	min4cost, min4edges, min4skipped = get_min_edges(normal_edge_cost, external_edge_cost, 4, cnt_skip)
	delta = cost10 - min4cost + 0.1 
	#print min4cost, costup, delta

	old_forbidden = forbidden

	forbidden = {}
	for i in range(N):
		forbidden[(i,i)] = True
		for j in range(i):
			if len(min4edges[i]) < 4 or len(min4edges[j]) < 4:
				continue # Not sure what to do in this case, be conservative
			mi = max(min4edges[i])[0]
			mj = max(min4edges[j])[0]
			if mi >= normal_edge_cost[i][j]:
				continue
			if mj >= normal_edge_cost[i][j]:
				continue
			# we are sure that normal_edge_cost[i][j] is more expensive then 4 best edges
			
			# min_increase (replace edge from min4edges by a fixed edge)
			min_increase = 2 * normal_edge_cost[i][j] - mi - mj
			
			# now, by fixing this edge, the skipped edges might have changed
			# Note that the fixed edge was not skipped in the first place (len(minedges4)==4)
			# so the only possibility is that mi and mj were skipped and currently cannot be
			# This is however good because new skipped length can be just smaller, i.e.
			# the real min_increase will be even bigger. Thus, this is a good estimate of the lower
			# bound
			if min_increase - 0.01 > delta:
				# by forcing this directed edge and removing maximum remaining edge
				# we still obtain worse solution
				forbidden[(i,j)] = True
				forbidden[(j,i)] = True
	cost_fast = solve_using_lp(normal_edge_cost, external_edge_cost, cnt_skip, forbidden)
	#print "normal:"
	CHECK_CONSISTENCY = False
	if CHECK_CONSISTENCY:
		cost_full = solve_using_lp(normal_edge_cost, external_edge_cost, cnt_skip, {})
		#print "check", cost_fast, cost_full
		assert abs(cost_fast - cost_full) < 0.0001
	
	#print "LB", min4cost, "UB", cost10, "fast", cost1, "norm", cost2
	#assert cost_fast <= cost10 + 0.0001
	if abs(cost10 - cost_fast) > 0.0001:
		print "10EDGES resulted in", cost10, "vs", cost_fast
	return cost_fast

def solve_using_lp(normal_edge_cost, external_edge_cost, cnt_skip, forbidden):
	start_time = time.time()
	external_selected = []
	for i in range(N):
		external_selected.append(Numberjack.VarArray(4, 0, 1))

	normal_selected = []
	for i in range(N):
		normal_selected.append(
			[ Numberjack.Variable(0, 1) if (i,j) not in forbidden else None 
				for j in range(i) ])
	vertex_skip = Numberjack.VarArray(N, 0, cnt_skip)
	
	model = Numberjack.Model()

	# ideme vytvorit constraint pre vrchol i
	for i in range(N):
		s = []
		s += external_selected[i]

		# edges j < i
		s += [ normal_selected[i][j] for j in range(i) if (i,j) not in forbidden]
		
		# edges j > i
		s += [ normal_selected[j][i] for j in range(i+1, N) if (i,j) not in forbidden]
		
		# we select specified number of outgoing edges
		model.add(Numberjack.Sum(s) == 4 - vertex_skip[i]);

	# limit number of edges to skip
	model.add(Numberjack.Sum(vertex_skip) <= cnt_skip)

	a = []
	b = []
	for i in range(N):
		for j in range(4):
			a.append(external_selected[i][j])
			b.append(external_edge_cost[i][j])
		for j in range(i):
			if (i,j) in forbidden:
				continue
			a.append(normal_selected[i][j])
			b.append(2 * normal_edge_cost[i][j])

	objective = Numberjack.Variable(-1e20, 1e20);
	model.add(Numberjack.Minimise(objective))
	model.add([objective == Numberjack.Sum(a, b)])
	solver = model.load('SCIP')
	solver.setVerbosity(0)
	res = solver.solve()
	assert res == True

	cost = objective.get_value()
	end_time = time.time()
	print "time", end_time - start_time
	solver.delete()
	return cost

cost = []
for cnt_skip in range(5):
	res = solve(normal_edge_cost, external_edge_cost, cnt_skip)
	cost.append(res);
fout = open(fout_filename, "w")
print >>fout, cost[0], cost[1], cost[2], cost[3], cost[4]
fout.close()
