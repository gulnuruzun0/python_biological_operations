'''
Ribosome Function with aminoacid dictionary
'''

aa = open("codon-aa-list.txt", "r")
list_aa = aa.readlines()[1:]

aa_list = []; counter = 0
for i in list_aa:
    aa_list.append(i.split('\t'))
    aa_list[counter][1] = aa_list[counter][1].replace("\n", "")
    counter += 1

aa.close()
aa_dict = {aa_list[i][0]: aa_list[i][1] for i in range(0, len(aa_list))} # dictionary comprehension
# print(aa_dict)

nucleotides = ["a", "c", "g", "u"]
seq = open("seq.txt", "r")
raw_seq = seq.read() # this method reads the whole file as a string not as a list of lines
raw_seq = [i.replace("\n", "") for i in raw_seq]
raw_seq = [i.replace(" ", "") for i in raw_seq]
raw_seq = [i.replace("\t", "") for i in raw_seq]
raw_seq = "".join([i for i in raw_seq if i in nucleotides]) 
raw_seq = raw_seq.upper()
# print(raw_seq)
# print(type(raw_seq))

seq.close()

def RNA_CDC(rna, aadict = aa_dict):
    start_index = ""
    gene = []
    for i in range(0, len(rna)):
        if rna[i:i+3] == "AUG":
            start_index = f"başlama indexi = {i+1}"
            gene.append(aadict[rna[i:i+3]])
            gene.append("-")
            for j in range(i+3, len(rna), 3):
                if rna[j:j+3] == "UAA" or rna[j:j+3] == "UAG" or rna[j:j+3] == "UGA":
                    gene.append(aadict[rna[j:j+3]])
                    return "".join(gene), start_index
                else:
                    gene.append(aadict[rna[j:j+3]])
                    gene.append("-")
            

    gene = "".join(gene)
    return gene, "there is no start codon"


# to readibility of the code
with open("AA_string.txt", "w") as gene: # "with" helps to close the file automatically 
    result = RNA_CDC(raw_seq)[0]
    for i in range(0, len(result), 60):
        if i + 60 < len(result):
            gene.write(result[i:i+60] + "\n")
        else:
            gene.write(result[i:])

# print(RNA_CDC(raw_seq)[0])
# print(RNA_CDC(raw_seq)[1])
