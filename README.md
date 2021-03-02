# test_scripts_PaulineFourgoux

## createmultifasta.py
It creates a multi fasta file from the lists of sequence files and blastp results given as arguments. It uses several arugments:
  * -i : input file : file containing the names of sequence and result files : seqNres.txt
  * -o : output file : name of output fasta file

Example command line call :
```
python3 createmultifasta.py -i <seqNres.txt> -o <outputname.fasta>
```

## simdimcount.py
It calculates the percentages of simmilarities and differences of an alignment of sequences. It uses several arugments:
  * -i : input file : result file from clustal omega : clustalo-seq_1-seq_2.txt
  * -f : fasta file which was the input for clustal omega : virus_seq.fasta

Example command line call :
```
python3 simdifcount.py -i <clustalo-res.txt> -f <sequences.fasta>
```
