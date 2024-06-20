'''
Alper Tanrıkulu, Gülnur Uzun

Q2:
(30 points) Missense Mutation Function
'''

from groupNo26_HW2_Q1 import RNA_CDC, aa_dict

def Mutation(seq:str, mutation:str, codons:dict = aa_dict): #obligatory parameters were defined along with the requested data types
    nucleotides = ["a", "c", "g", "u"]
    seq = [i.replace("\n", "") for i in seq]
    seq = [i.replace(" ", "") for i in seq]
    seq = [i.replace("\t", "") for i in seq]
    seq = "".join([i for i in seq if i in nucleotides]) 
    seq = seq.upper() #with string operation, the string in fasta format was converted into a straight RNA sequence written in large character

    position = ''.join(filter(str.isdigit, mutation)) 
    
    #mutation information has been translated into the list to get the mutated sequence
    # mut = [mutation[0], mutation[2:5], mutation[5], mutation[7]] 
    mut = [mutation[0], position, mutation[-3], mutation[-1]] 
    
    #in order to edit on it, the reference array was converted to a list and assigned to a new variable
    mutated_seq = list(seq) 
    mutated_seq[int(mut[1])-1] = str(mut[3]) #the given sequence has been converted to mutated
    mutated_seq = "".join(mutated_seq) #the mutated list has been converted to a string.

    mutated = RNA_CDC(mutated_seq)[0] #with the function we defined earlier, the reference sequence was converted to the aa sequence.
    not_mutated = RNA_CDC(seq)[0] #the mutated sequence was converted to the aa sequence
    print(not_mutated)
    for i in range(0,max(len(mutated), len(not_mutated)), 4): #to compare the strings from the function, the loop is rotated through the longest array.
        #there is a "-" sign between amino acids in the string, which is the output of the function, for the purpose of convenient decodability, and the for loop was created taking this into account
        
        if mutated[i:i+3] != not_mutated[i:i+3]: #if different amino acids are seen in the same position, that is, if there is a missense mutation;
            return f"p.{not_mutated[i:i+3]}{int(i/4)+1}{mutated[i:i+3]}" #returning mutation information in the desired format
    return f"The mutation is not missense mutation!"



#---------------------------------------------------------------------------------------------------------------------------------------
# You can control the parameters here.
mutat1 = "r.69A>C"
fasta_seq = open("seq.txt", "r")
seq = fasta_seq.read()

print(Mutation(seq, mutat1))



