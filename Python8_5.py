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
                seqs[identificador]["rev_comp_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC")) 
                # adiciona a sequencia reversa complementar ao dicionario
            id = re.search(r'>(\S+)(\s.+?)', l) # busca o identificador
            identificador = id.group(1)
            seqs[identificador] = {"seq": "", "frame_+1":[], "frame_+2":[], "frame_+3":[], "rev_comp_seq": "", "frame_-1":[], "frame_-2":[], "frame_-3":[]} # cria um dicionario para o identificador
        else:
            seqs[identificador]["seq"] += l.upper() # adiciona a sequencia ao dicionario
    seqs[identificador]["rev_comp_seq"] = seqs[identificador]["seq"][::-1].translate(str.maketrans("ATCG", "TAGC"))
    
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
                
            
            ##Reverse

            frame = "frame_-" + str(i+1) 

            for match in re.finditer(r"(.{3})", seqs[id]["rev_comp_seq"][i:]): # busca todas as bases em grupos de 3
                seqs[id][frame].append(match.group(1)) # adiciona ao dicionario
                
            if len(seqs[id][frame][-1]) != 3: # o último elemento pode ter menos de 3 bases, então removemos se for o caso
                seqs[id][frame][-1].pop() # remove o ultimo elemento
            
            # salvando os resultados em formato fasta, na ordem de id e frames (1 a 3 e -1 a -3)
            out += f"{id}-frame-{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
            out_rev += f"{id}-frame--{i+1}-codons\n{' '.join(seqs[id][frame])}\n"
        
        # escrevendo os resultados
        f.write(out)
        f.write(out_rev)

translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

with open("Python_08.translated.aa", "w") as f:
    # traduzindo os frames para proteina
    for id in seqs:
        for i in range(3): 
            frame = "frame_+" + str(i+1)
            protein_frame = "protein_frame_+" + str(i+1)
            seqs[id][protein_frame] = "" #criando uma chave para armazenar a proteina em cada frame

            for codon in seqs[id][frame]: #para cada codon no frame
                if not codon in translation_table: #se o codon não estiver na tabela de tradução, substituir por _
                    print(f"Warning: codon {codon} not found in translation table, using _ as placeholder")
                    seqs[id][protein_frame] = "_"
                else:
                    seqs[id][protein_frame] += translation_table[codon]

            ## Reverse, mesma coisa que o anterior, mas para o frame negativo

            frame = "frame_-" + str(i+1)
            protein_frame = "protein_frame_-" + str(i+1)
            seqs[id][protein_frame] = ""

            for codon in seqs[id][frame]:
                if not codon in translation_table:
                    print(f"Warning: codon {codon} not found in translation table, using _ as placeholder")
                    seqs[id][protein_frame] = "_"
                else:
                    seqs[id][protein_frame] += translation_table[codon]

        #salvando os resultados em formato fasta, na ordem de id e frames (1 a 3 e -1 a -3)
        f.write("".join([f">{id}-frame-{i}-protein\n{seqs[id]['protein_frame_+'+str(i)]}\n" for i in range(1,4)]))
        f.write("".join([f">{id}-frame--{i}-protein\n{seqs[id]['protein_frame_-'+str(i)]}\n" for i in range(1,4)]))

