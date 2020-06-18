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

def SequenceToString(fasta_sequence):
    #returns sequence as str

    final_str = str(fasta_sequence)
    return final_str

def FixOrganismName(old_organism_name):
    # removes any parenthesis with common name

    str2 = re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", old_organism_name)

    final_str = str2.replace("()", "")

    return final_str


numberRecords("arc_sequences_06172020.gp")

file = open("output_repeated_sequences.csv", "w")
file.write("#" + "," + "Source" + "," + "# Amino Acids" +  "," + "Sequence" + "," + "Comments" + "\n")

## Creates a list of all sequences in list format
ListOfSequences = []

for seq_record in SeqIO.parse("arc_sequences_06172020.gp", "gb"):
    SeqToAppend = str(seq_record.seq)
    ListOfSequences.append(SeqToAppend)

#################################################

def CheckWholeList(ListRecords, ListOfSequences):
    #assigns group, write with sequence to a file, (in progress) remove duplicate sequences or unknown XXXX

    counter = 0
    counterRecs = 0
    old_sequence = "no sequence yet"

    for seq_record in SeqIO.parse(ListRecords, "gb"):  #for every record in the list

        # setting up initial vatiables
        new_sequence_name = seq_record.annotations["source"]
        new_sequence_length = len(seq_record)
        new_sequence = str(seq_record.seq)


        counterRecs = counterRecs + 1  #counter of records
        recNumber = str(counterRecs)
        printSequence = SequenceToString(seq_record.seq)
        printSciName = FixOrganismName(seq_record.annotations["source"])
        comment = " "

        if (new_sequence == old_sequence):
            comment = "Duplicates previous entry  "

        ct = ListOfSequences.count(printSequence)

        if ct > 1:
            comment3 = "Sequence appears total of: " + str(ct) + " times in the dataset"
            file.write(recNumber + "," + printSciName + "," + str(new_sequence_length) + "," + printSequence + "," + comment + "," + comment3 + "\n")

        old_sequence = new_sequence

    print("Number of unclassified species:", counter)
    print("Number of records written to file: ", counterRecs)

    file.close()

CheckWholeList("arc_sequences_06172020.gp", ListOfSequences)