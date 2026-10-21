import networkx as nx,json,time
j=json.load(open('data/enumerations/vinberg_x3_h19.json'));G=j['wall_gram'];g=nx.Graph();g.add_nodes_from((i,{'norm':G[i][i]})for i in range(36));g.add_edges_from((i,k,{'pair':G[i][k]})for i in range(36)for k in range(i)if G[i][k])
gm=nx.algorithms.isomorphism.GraphMatcher(g,g,node_match=lambda a,b:a==b,edge_match=lambda a,b:a==b)
a=[];t=time.time()
for p in gm.isomorphisms_iter():
 a.append([p[i]for i in range(36)])
 if len(a)%1000==0:print(len(a),time.time()-t,flush=True)
 if len(a)>=20000:break
json.dump({'automorphisms':a,'limit':20000,'exhausted':len(a)<20000},open('data/enumerations/fresh_h20_coupled_completion_20261013_automorphisms.json','w'),indent=2);print(len(a),time.time()-t,flush=True)
