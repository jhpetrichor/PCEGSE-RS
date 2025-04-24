import networkx as nx


def read_essential_protien() -> set:
    lines = open('essential proteins.ref', 'r').readlines()[:]
    ess = set()
    
    for line in lines:
        line = line.splitlines()[0]
        ess.add(line)

    return ess

def read_ppi_protein(ppi: str):
    lines= open(ppi, 'r').readlines()
    protein = set()
    for line in lines:
        print(line)
        line = line.split()[:]
        protein.add(line[0])
        protein.add(line[1])
    return protein

def get_label(essential: set, protein: set):
    label = {}
    for p in protein:
        if p in essential:
            label[p] = 1
        else:
            label[p] = 0
    return label

ess = read_essential_protien()
data = ['BioGRID_Y2H', 'IntAct_AP-MS', 'IntAct', 'BioGRID_Y2H']
for d in data:
    ppi = d + "/" + d + '.txt'
    protein = read_ppi_protein(ppi)
    label = get_label(ess, protein)
    path = d + '/' + d + '_label.txt'
    file = open(path, 'w')
    for a in label.items():
        file.write(a[0] + '\t' + str(a[1]) + '\n')

for d in data:
    G = nx.read_edgelist(d + "/" + d + '.txt')
    cliques = list(nx.find_cliques(G))

    # 将团写入文件
    output_file = d + '/' + d + '_clique.txt'
    with open(output_file, "w") as f:
        for clique in cliques:
            for p in clique:
                f.write(f"{p}\t")
            f.write(f"\n")

    print(f"所有团已写入 {output_file}")