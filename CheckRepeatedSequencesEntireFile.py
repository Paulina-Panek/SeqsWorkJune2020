# Paulina Panek
# June 2020
# Script goes one by one sequence and compares it against others to see if there is another identical sequence elsewhere

from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO
import re

def numberRecords(ListRecords):
    #Function prints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records in initial file " % len(records))

def CheckIfDuplicate(first_sequence, second_sequence):
    # returns 0 (same sequences), 1 (not same sequences, or 3 (something went wrong, function didn't work)

    return_value = 3

    # if same species AND length of sequence is the same, check if the sequence is the same
    if first_sequence == second_sequence:
        return_value = 0  #same sequences
    else:
        return_value = 1

        return(return_value)

def unknown_aas(sequence):
    #returns number of unknown amino acids in sequence
    X_in_sequence = 0

    if 'X' in sequence:
        X_in_sequence = X_in_sequence + 1
    return X_in_sequence

def SequenceToString(fasta_sequence):
    #returns sequence as str

    final_str = str(fasta_sequence)
    return final_str

def FixProteinName(old_protein_name):
    # removes , from protein name (which confuses Excel) and removes any parenthesis
    str1 = old_protein_name.replace(",", ";")

    str2 = re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", str1)

    final_str = str2.replace("[]", "")

    return final_str


numberRecords("arc_sequences_06172020.gp")

file = open("output_repeated_sequences.txt", "w")
file.write("#" + "\t" + "Protein Name" + "\t" + "# Amino Acids" + "\t" + "Accession Number" + "\t" + "Source Organism" +  "\t" + "Sequence" + "\t" + "Comment" + "\n")


def MakeExcel(ListRecords):
    #assigns group, write with sequence to a file, (in progress) remove duplicate sequences or unknown XXXX

    counter = 0
    counterRecs = 0
    old_sequence = "no sequence yet"

    for seq_record in SeqIO.parse(ListRecords, "gb"):  #for every record in the list

        # setting up initial vatiables
        new_sequence_name = seq_record.annotations["source"]
        new_sequence_length = len(seq_record)
        new_sequence = str(seq_record.seq)
        Number_of_X = unknown_aas(new_sequence)

        counterRecs = counterRecs + 1  #counter of records
        recNumber = str(counterRecs)
        printname = FixProteinName(seq_record.description)
        printSequence = SequenceToString(seq_record.seq)
        comment = " "

        if (CheckIfDuplicate(new_sequence, old_sequence) == 0):
            comment = ("*REMOVED* Duplicate  " )
        elif (Number_of_X != 0):
            comment = ("*REMOVED* Unknown AAs  ")

        file.write(recNumber + "\t" + printname + "\t" + str(new_sequence_length)  + "\t" + seq_record.id + "\t" + seq_record.annotations["source"] + "\t" + printSequence + "\t" + comment +"\n")

        old_sequence = new_sequence

    print("Number of unclassified species:", counter)
    print("Number of records written to file: ", counterRecs)

    file.close()

MakeExcel("arc_sequences_06172020.gp")


#for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
#print(seq_record.description) #protein name [organism]
#print(seq_record.seq) # sequence
#print(seq_record.annotations["source"]) #name (common name)
#print(seq_record.annotations["taxonomy"][0])
#print(seq_record.annotations)
#print(len(seq_record)) #length of sequence
