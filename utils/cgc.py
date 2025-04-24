import networkx as nx

def read_ppi(file: str, weight=False):
    edges = []
    nodes = set()
    lines = open(file, 'r').readlines()
    for line in lines:
        edge = line.strip().split('\t')[:]
        # print(edge)
        nodes.add(edge[0])
        nodes.add(edge[1])
        if weight:
            edges.append((edge[0], edge[1], float(edge[2])))
        else:
            edges.append((edge[0], edge[1]))
    return edges, nodes


import networkx as nx


def split(queue: list[nx.Graph], splitted: list[nx.Graph]):
    G = queue.pop()
    n = G.number_of_nodes()

    if n < 20:
        if nx.is_connected(G):
            splitted.append(G)
        else:
            for comp in nx.connected_components(G):
                subgraph = G.subgraph(comp).copy()
                splitted.append(subgraph)
    else:
        if not nx.is_connected(G):
            for comp in nx.connected_components(G):
                subgraph = G.subgraph(comp).copy()
                queue.append(subgraph)
        else:
            # 对图进行 Kernighan–Lin 二分
            try:
                part1, part2 = nx.algorithms.community.kernighan_lin_bisection(G)
                G1 = G.subgraph(part1).copy()
                G2 = G.subgraph(part2).copy()
                queue.append(G1)
                queue.append(G2)
            except Exception as e:
                # 万一切不了（比如图太小或结构不支持），也保底放进结果
                splitted.append(G)


edges, _ = read_ppi('../data/collins/collins.txt', False)

G = nx.from_edgelist(edges)

splitted = []
queue = [G]

while len(queue) > 0:
    split(queue, splitted)

for g in splitted:
    if len(g.nodes()) < 3:
        continue
    for i in g.nodes():
        print(i, end=' ')
    print()
