#!usr/bin/python3

#---------------------------------------------------#
# PACKAGES CALL
#---------------------------------------------------#
import math
import sys

#---------------------------------------------------#
# MAIN
#---------------------------------------------------#
def main():
    # Retrieve arguments
    input_name, output_name = getArgs()

    # Read input file
    input = open(input_name, 'r')
    inputlines = input.readlines()
    input.close()

    # Open output file
    output = open(output_name, 'w')

    # For each line of input file
    for i in range(len(inputlines)):
        line = inputlines[i]
        line = line.rstrip() # remove \n character
        line = line.split("\t")
        header = readBlastResults(line[1])
        fasta = txttofasta(line[0])
        output.write(header + "\n")
        for i in range(len(fasta)):
            output.write(fasta[i] + "\n")

    print("Multi-fasta file created.")



#---------------------------------------------------#
# FUNCTIONS
#---------------------------------------------------#
def usage():
    print("""
    createmultifasta.py creates a multi fasta file from the lists of sequence
    files and blastp results given as arguments.

    exemple of command line:
    =======================

    python3 createmultifasta.py -i seqNres.txt -o outputname.fasta

    obligatory:
    ==========

    -i -> file containning the names of sequence and result files

    -o -> name of output fasta file

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
        output_file = sys.argv[sys.argv.index("-o")+1]
        print("output file : ", output_file)
    except:
        usage()
        print("ERROR : please enter the output file")
        sys.exit()

    return (input_file, output_file)


def txttofasta(txtfile):
    '''
        txttofasta a function that takes a file name of a character string as
        argument and returns the list of lines for a fasta file comprised of
        80 character or less.
    '''
    # Create empty list of fasta lines
    fasta = []
    # Open sequence file
    filetoread = open(txtfile)
    # Retrieve the unique line with the sequence
    str_seq = filetoread.readline()
    # Close file
    filetoread.close()
    str_seq = str_seq.rstrip() # remove \n character
    lenseq = len(str_seq)
    # Retrieve the number of full fasta lines (80 characters) possible from the sequence
    nb80full = math.floor(lenseq/80)

    for i in range(nb80full):
        # Create a fasta line of 80 character
        fasta.append(str_seq[(i*80):(i*80+80)])

    # Add the last fasta line with the remaining characters
    fasta.append(str_seq[(nb80full*80):lenseq])

    return(fasta)

def readBlastResults(blastfilename):
    '''
        readBlastResults is a function that takes as argument the file name of
        blast results and take the first line of results which is the one with
        the best score and returns the header for a fasta file.
    '''
    # Open file and read file
    resfile = open(blastfilename)
    res = resfile.readlines()
    resfile.close()
    # Read first result line with the best score (line 2)
    line2 = res[1]
    line2 = line2.rstrip() # remove \n character
    line2 = line2.split(",") # split the line by the comma for column separation
    access = line2[9].split('\"')[2]
    fastaheader = ">" + access + " " + line2[0]
    return fastaheader

#---------------------------------------------------#

if __name__ == "__main__":
    main()
