import re
import sys
import os

seqs = {} # dicionario para armazenar os identificadores como chave e as bases como valor

if len(sys.argv) < 2:
    print("No file")
    sys.exit(1)

file = sys.argv[1]

if not os.path.exists(file):
    print("File not found")
    sys.exit(1)

with open(file) as f:
    for l in f:
        l = l.rstrip().upper() # remove \n e coloca em maiuscula, para evitar erros
        if l.startswith(">"): # testa se é cabeçalho
            id = re.search(r'>(\S+)(\s.+?)', l) # busca o identificador
            identificador = id.group(1) # pega o identificador
            seqs[identificador] = {"A":0, "T":0, "C":0, "G":0} # cria um dicionario para o identificador
        else:
            seqs[identificador]["A"] += l.count("A") # conta as bases e adiciona ao dicionario
            seqs[identificador]["T"] += l.count("T")
            seqs[identificador]["C"] += l.count("C")
            seqs[identificador]["G"] += l.count("G")
    
for id in seqs:
    print(f"{id}\tA_{seqs[id]['A']}\tT_{seqs[id]['T']}\tC_{seqs[id]['C']}\tG_{seqs[id]['G']}")
