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
    identificador = ""
    for l in f:
        l = l.rstrip()
        if l.startswith(">"):  # testa se é cabeçalho
            if identificador:
                seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
            id = re.search(r'>(\S+)(\s.+?)', l) # busca o identificador
            identificador = id.group(1)
            seqs[identificador] = {"seq": "", "frame_+1":[], "frame_+2":[], "frame_+3":[], "reverse_seq": "", "frame_-1":[], "frame_-2":[], "frame_-3":[]} # cria um dicionario para o identificador
        else:
            seqs[identificador]["seq"] += l.upper() # adiciona a sequencia ao dicionario
    seqs[identificador]["reverse_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
    
with open("Python_08.codons-6frames.nt", "w") as f:

    for id in seqs: # para cada identificador
        out = ""
        out_rev = ""

        for i in range(3): # para cada frame, no caso só o +1
        
            frame = "frame_+" + str(i+1) 

            for match in re.finditer(r"(.{3})", seqs[id]["seq"][i:]): # busca todas as bases em grupos de 3
                seqs[id][frame].append(match.group(1)) # adiciona ao dicionario
                
            if len(seqs[id][frame][-1]) != 3: # o último elemento pode ter menos de 3 bases, então removemos se for o caso
                seqs[id][frame][-1].pop() # remove o ultimo elemento
            out += f"{id}-frame-{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
            
            ##Reverse

            frame = "frame_-" + str(i+1) 

            for match in re.finditer(r"(.{3})", seqs[id]["reverse_seq"][i:]): # busca todas as bases em grupos de 3
                seqs[id][frame].append(match.group(1)) # adiciona ao dicionario
                
            if len(seqs[id][frame][-1]) != 3: # o último elemento pode ter menos de 3 bases, então removemos se for o caso
                seqs[id][frame][-1].pop() # remove o ultimo elemento
            
            out_rev += f"{id}-frame--{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
        
        f.write(out)
        f.write(out_rev)

