import networkx as nx

# 创建一个图
G = nx.Graph()

# 添加边
G.add_edges_from([
    (1, 2),
    (1, 3),
    (2, 3),
    (2, 4),
    (3, 4),
    (4, 5)
])

# 找到所有团
cliques = list(nx.find_cliques(G))

# 打印所有团
print("图中的所有团：")
for clique in cliques:
    print(clique)