#!usr/bin/python3

#---------------------------------------------------#
# PACKAGES CALL
#---------------------------------------------------#
import sys

#---------------------------------------------------#
# MAIN
#---------------------------------------------------#
def main():
    # Retrieve arguments
    input_name, fasta_name = getArgs()

    # Retrieve accession ID of sequences aligned
    seqs = retrieveseqID(fasta_name)

    # Counts
    conservedRes, strongConserv, weakConserv, noConserv, total = simdifcount(input_name, seqs)
    print("Alignment similarities and differences :")
    print("Percentage of conserved residues : ", round(conservedRes/total*100,2), "%")
    print("Percentage of strongly conserved residue properties : ", round(strongConserv/total*100,2), "%")
    print("Percentage of weakly conserved residue properties : ", round(weakConserv/total*100,2), "%")
    print("Percentage of residues not conserved : ", round(noConserv/total*100,2), "%")
    print("Percentage of overall conservation : ", round((conservedRes+strongConserv+weakConserv)/total*100,2), "%")

#---------------------------------------------------#
# FUNCTIONS
#---------------------------------------------------#
def usage():
    print("""
    simdifcount.py calculates the percentages of simmilarities and differences
    of an alignment of sequences.

    exemple of command line:
    =======================

    python3 simdifcount.py -i clustalo-res.txt

    obligatory:
    ==========

    -i -> result file from clustal omega

    -f -> fasta file which was the input for clustal omega

    """)

def getArgs():
    try:
        input_file = sys.argv[sys.argv.index("-i")+1]
        print("file to treat : ", input_file)
    except:
        usage()
        print("ERROR : please enter the input file")
        sys.exit()
    try:
        fasta_file = sys.argv[sys.argv.index("-f")+1]
        print("fasta to treat : ", fasta_file)
    except:
        usage()
        print("ERROR : please enter the fasta file")
        sys.exit()

    return (input_file, fasta_file)


def retrieveseqID(fasta):
    fastafile = open(fasta, 'r')
    fastalines = fastafile.readlines()
    fastafile.close()

    # Lines starting with '>'
    res = []
    for i in range(len(fastalines)):
        if fastalines[i][0] == '>':
            ID = fastalines[i].split(" ")[0]
            ID = ID[1:] # remove the >
            res.append(ID)

    return res




def simdifcount(clustalres, seq):
    '''
        Results from clustal omega meaning
        * -> indicates positions which have a single, fully conserved residue.
        : -> indicates conservation between groups of strongly similar properties as below - roughly equivalent to scoring > 0.5 in the Gonnet PAM 250 matrix.
        . -> indicates conservation between groups of weakly similar properties as below - roughly equivalent to scoring =< 0.5 and > 0 in the Gonnet PAM 250 matrix.
          -> indicates a gap or a non conservation.
    '''

    # Create a dictionnary for each sequence ID to reform the sequence
    dico_seq = {}

    for access in seq:
        dico_seq[access] = ''
    dico_seq["alignment"] = ''

    clustalfile = open(clustalres, 'r')
    clustallines = clustalfile.readlines()
    clustalfile.close()

    # Remove header from the Results
    clustallines = clustallines[3:]

    # calculate the number of caracter befor begining of sequences
    start = len(clustallines[0].split("\t")[0]) - 60

    for i in range(len(clustallines)):
        if (i+1)%4 != 0: # Every 4 lines there is an empty line, just with a \n character
            line = clustallines[i].strip("\n") # remov \n character of the string
            ID = line.split(" ")[0]
            line = line.split("\t")[0]
            sequence = line[start :] # retrieving only the 60 characters of the alignment as clustal omega only print 60 characters per line
            if len(ID)>0:
                dico_seq[ID] += sequence
            else:
                dico_seq["alignment"] += sequence

    conservedRes = dico_seq["alignment"].count('*')
    strongConserv = dico_seq["alignment"].count(':')
    weakConserv = dico_seq["alignment"].count('.')
    noConserv = dico_seq["alignment"].count(' ')
    total = conservedRes + strongConserv + weakConserv + noConserv

    return conservedRes, strongConserv, weakConserv, noConserv, total

#---------------------------------------------------#

if __name__ == "__main__":
    main()
