# Paulina Panek
# June 2020
# Script parsing result GenPept file (.gp) to get a csv file with all results

from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO
import re

def numberRecords(ListRecords):
    #Function prints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records in initial file " % len(records))

def CheckIfDuplicate(first_sequence_name, second_sequence_name, first_sequence, second_sequence):
    # returns 0 (same sequences), 1 (not same sequences, or 3 (something went wrong, function didn't work

    return_value = 3

    # if same species AND length of sequence is the same, check if the sequence is the same
    if (first_sequence_name == second_sequence_name):
        if first_sequence == second_sequence:
            return_value = 0  #same sequences
        else:
            return_value = 1
    else:
        return_value = 1
    return(return_value)

def unknown_aas(sequence):
    #returns number of unknown amino acids in sequence
    X_in_sequence = 0

    if 'X' in sequence:
        X_in_sequence = X_in_sequence + 1
    return X_in_sequence

def TaxonomyToString(Taxonomy, separator):
    #Function converts dict in the taxonomy field to one string

    final_str = separator.join(Taxonomy)
    return final_str

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

file = open("output_masterfile.csv", "w")
file.write("#" + "," + "Protein Name" + "," + "# Amino Acids" + "," + "Accession Number" + "," + "Source Organism" + "," + "Full Taxonomy" + "," + "Assignment" + "," + "Sequence" + "," + "Comment" + "\n")


def MakeExcel(ListRecords):
    #assigns group, write with sequence to a file, (in progress) remove duplicate sequences or unknown XXXX

    counter = 0
    counterRecs = 0
    duplicates = 0
    old_sequence_name = "empty"
    old_sequence = "no sequence yet"
    new_sequence_name = "empty2"
    sequence_title = "error! check what happened here"

    for seq_record in SeqIO.parse(ListRecords, "gb"):  #for every record in the list

        # setting up initial vatiables
        new_sequence_name = seq_record.annotations["source"]
        new_sequence_length = len(seq_record)
        new_sequence = str(seq_record.seq)
        assignment = "UNASSIGNED FIX ME"
        Number_of_X = unknown_aas(new_sequence)
        prot_name = seq_record.description

        #if (CheckIfDuplicate(new_sequence_name, old_sequence_name, new_sequence, old_sequence) == 1) and (Number_of_X == 0) and RemoveLike(prot_name) == 0:  # if not the same and no unknown aas (X), continue


#Classification block begins~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (seq_record.annotations["taxonomy"][2] == "Ecdysozoa") or (seq_record.annotations["taxonomy"][4] == "Insecta"):  #classify as invertebrate
            assignment = "Invertebrate"

        elif seq_record.annotations["taxonomy"][6] == "Amphibia":  # classify as amphibia
            assignment = "Amphibian"

        elif (seq_record.annotations["taxonomy"][6] == "Actinopterygii") or (seq_record.annotations["taxonomy"][6] == "Coelacanthiformes"):  # classify as fish
            assignment = "Fish"

        elif seq_record.annotations["taxonomy"][6] == "Archelosauria":  # classify as reptile or bird
            if seq_record.annotations["taxonomy"][11] == "Coelurosauria" or seq_record.annotations["taxonomy"][11] == "Aves": #bird
                assignment = "Bird"
            else:
                assignment = "Reptile"

        elif seq_record.annotations["taxonomy"][6] == "Archosauria":  # classify as bird
            if seq_record.annotations["taxonomy"][11] == "Aves":  # bird
                assignment = "Bird"
            else:
                counter = counter + 1

        elif seq_record.annotations["taxonomy"][6] == "Lepidosauria" or  seq_record.annotations["taxonomy"][6] == "Testudines + Archosauria group":
            assignment = "Reptile"

        elif seq_record.annotations["taxonomy"][6] == "Mammalia":
            if seq_record.annotations["taxonomy"][9] == "Primates":
                assignment = "Primate"
            else:
                assignment = "Mammal"
        else:
            assignment = "UNCLASSIFIED FIX ME"
            #counter = counter + 1
#end of classification block~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        counterRecs = counterRecs + 1  #counter of records
        recNumber = str(counterRecs)
        printname = FixProteinName(seq_record.description)
        printTaxonomy = TaxonomyToString(seq_record.annotations["taxonomy"], ";")
        printSequence = SequenceToString(seq_record.seq)
        comment = " "

        if (CheckIfDuplicate(new_sequence_name, old_sequence_name, new_sequence, old_sequence) == 0):
            comment = ("*REMOVED* Duplicate  " )
        elif (Number_of_X != 0):
            comment = ("*REMOVED* Unknown AAs  ")

        file.write(recNumber + "," + printname + "," + str(new_sequence_length)  + "," + seq_record.id + "," + seq_record.annotations["source"] + "," + printTaxonomy + "," + assignment + "," + printSequence + "," + comment +"\n")

        old_sequence_length = new_sequence_length
        old_sequence_name = new_sequence_name
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
