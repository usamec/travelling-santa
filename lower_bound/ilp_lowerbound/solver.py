#!/usr/bin/python
import sys
sys.path.append("/home/ppershing/python/numberjack/local_lib");
sys.setrecursionlimit(4000)
import Numberjack

fin = open("solver.in", "r")
tmp = fin.readline()
tmp = fin.readline();
N = int(tmp.strip())
print "Problem dimension: ", N
external_edges = []
for i in range(N):
	tmp = fin.readline()
	tmp = tmp.strip().split(' ')
	assert len(tmp) == 4
	external_edges.append([])
	for j in range(4):
		external_edges[i].append(float(tmp[j]))


normal_edges = []
for i in range(N):
	tmp = fin.readline()
	tmp = tmp.strip().split(' ')
	assert(len(tmp) == N)
	normal_edges.append([])
	for j in range(N):
		normal_edges[i].append(float(tmp[j]))

fin.close()

for i in range(N):
	for j in range(N):
		assert(normal_edges[i][j] == normal_edges[j][i])

cost = []
for edges_to_skip in range(5):
	external_selected = []
	external_skip = []
	for i in range(N):
		external_selected.append(Numberjack.VarArray(4, 0, 1))
		external_skip.append(Numberjack.VarArray(4, 0, 1))

	normal_selected = [None,]
	normal_skip = [None,]
	for i in range(1, N):
		normal_selected.append(Numberjack.VarArray(i, 0, 1))
		normal_skip.append(Numberjack.VarArray(i, 0, 2)) # can skip in both directions

	model = Numberjack.Model()

	for i in range(N):
		# ideme vytvorit constraint pre vrchol i
		s = 0
		# external edges
		for j in range(4):
			s += external_selected[i][j]
		# edges j < i
		for j in range(i):
			s += normal_selected[i][j]
		# edges j > i
		for j in range(i+1, N):
			s += normal_selected[j][i]
		model.add(s >= 4); # we select at least 4 edges from each vertex

	# limit number of edges to skip
	skipped = 0;
	for i in range(N):
		skipped += Numberjack.Sum(external_skip[i])
		if i:
			skipped += Numberjack.Sum(normal_skip[i])
	model.add(skipped <= edges_to_skip)

	# can skip only if selected
	for i in range(N):
		for j in range(4):
			model.add(external_skip[i][j] <= external_selected[i][j]);
		for j in range(i):
			model.add(normal_skip[i][j] <= normal_selected[i][j] + normal_selected[i][j]);

	a = []
	b = []
	for i in range(N):
		for j in range(4):
			a.append(external_selected[i][j] - external_skip[i][j])
			b.append(external_edges[i][j])
		for j in range(i):
			a.append(normal_selected[i][j] +  normal_selected[i][j] - normal_skip[i][j])
			b.append(normal_edges[i][j])

	objective = Numberjack.Variable(-1e20, 1e20);
	model.add(Numberjack.Minimise(objective))
	model.add([objective == Numberjack.Sum(a, b)])
	solver = model.load('SCIP')
	solver.setVerbosity(1)
	res = solver.solve()
	assert res == True

	cost.append(objective.get_value())
fout = open("solver.out", "w")
print >>fout, cost[0], cost[1], cost[2], cost[3], cost[4]
fout.close()
