
---------------------------------------------------------------------------------
About the fstage
---------------------------------------------------------------------------------

The fstage software suite is a convenient tool for working with large amounts 
of biological data represented as nucleotide sequences. The initial purpose is 
to solve the problems of selecting consensus sequences for which primer design 
could be carried out in the future. For this, the program performs a wide range 
of functions: automation of routine launches of programs for multiple alignment 
(GramAlign, Clustalo, Muscle), grouping, filtering and clustering of nucleotide 
sequences; converting fasta alignments to consensus sequences; refinement of 
consensus sequences based on a set of short sequences.

-----------------------------------------------------------------------------------
Installation
-----------------------------------------------------------------------------------

The program is fully written in Python 3 and is supported on Windows and Linux. To 
run the program, you must first install the Python interpreter and the BioPython 
library.

Installation can be executed by two ways:

1. Cloning from a remote repository using the git version control system with the 
command:
                 git clone https://github.com/Maximato/fstage

2. Copying files in the working directory manually from the source 
                      https://github.com/Maximato/fstage

-----------------------------------------------------------------------------------
Script align.py
-----------------------------------------------------------------------------------

The program is designed to run multiple sequence alignment of one of the three 
programs (GramAlign, muscle or clustalo). These programs should be installed 
on your computer and added to path variable for correct work. If programs were 
not be added to 'path' system variable they can be added to definitions.py. 
The launch is as follows:

	python3 align.py -i sequences.fasta –p program_name -o out_file.fasta -m

-i, --input (sequences.fasta) - a set of sequences for which multiple alignment 
is necessary (or the name of the directory with the files with sequences 
stored in it);
-p, --program (program_name) - (optional parameter, default GramAlign) 
name of the program that performs alignment: GramAlign, muscle or clustalo
-o, --output (out_file.fasta) - name of the output file (or directory);
–m, --multiply - an optional parameter for automatically starting alignments 
of all files contained in the directory specified in the -i (--input) parameter.


-----------------------------------------------------------------------------------
Script alwork.py
-----------------------------------------------------------------------------------

This program is necessary for various manipulations with alignments in fasta format. 
Available commands: convert - converts aligned sequences to consensus, 
consensus - calculates a new consensus based on alignment of complete genomes and 
shorter sequences, unite - combines alignments into a single file containing 
consensus sequences.

Using convert:

	python3 alwork.py convert –i align.fasta –o out_dir –p prefix

-i, --input (align.fasta) - the result of alignment of sequences in fasta format;
-o, --output (out_dir) - the name of the directory for storing the calculation 
results;
-p, --prefix (prefix) - prefix file names in the results.

Using consensus:

	python3 alwork.py consensus –i align.fasta –s seqs.fasta –o output –f format

-i, --input (align.fasta) - the result of alignment of sequences in fasta format;
-s, --seqs (seqs.fasta) - file with short sequences in fasta format;
-o, --output (output) - name of the output file containing the new consensus;
-f, --format (format) - output file format (“html” or “fasta”).

Using unite:

	python3 alwork.py unite –i input –o output –fl False –ig True –il 0.9

-i, --input (input) - name of the directory containing files with alignments 
in fasta format;
-o, --output (output) - the name of the output file, with consensus sequences;
-fl, --flength - (optional parameter, default True) True or False, whether 
the ends of sequences consisting only of gaps (gaps) will be taken into account 
to calculate the occurrence of a nucleotide at a given alignment position;
-ig - (optional parameter, default False) True or False, ignoring gaps when 
drawing up a consensus sequence;
-il - (optional parameter, default 0.9) number, occurrence level above which 
gaps are ignored.

-----------------------------------------------------------------------------------
Script opseq.py
-----------------------------------------------------------------------------------

This program is designed to work with nucleotine sequences. With its help are 
available: grouping (group), filtering (filtr) and clustering (clust) sequences.

Using group:

	python3 opseq.py group –i sequences.fasta –o out_dir –minsog 100 –maxsog 300

-i, --input (sequences.fasta) - a file in fasta format containing the sequences to be 
studied;
-i, --output (out_dir) - the name of the directory in which the grouped sequences 
will be stored;
--minsog - non-negative integer, minimal genome size;
--maxsog - non-negative integer, the maximum size of the genome.

Using filtr:

	python3 opseq.py filtr –i sequences.fasta –o outfile –organism TBEV –mins 100 –maxs 300

-i, --input (sequences.fasta) - a file in fasta format containing the sequences 
to be studied;
-o, --output (outfile) - the name of the output file;
--organism (organism) - the full name of the organism for which filtering 
is performed;
--mins - non-negative integer, minimum sequence size;
--maxs - non-negative integer, the maximum size of the sequences.

Using clust:

	python3 opseq.py clust –i sequences.fasta –o out_dir –e 0.5 –s 2 --dm dist_matrix

-i, --input (sequences.fasta) - a file in fasta format containing the sequences 
to be studied;
-o, --output (out_dir) - the name of the directory in which clusters with 
sequences and visualization data will be stored;
-e, --eps - (optional parameter, default 0.5) non-negative number, maximum distance 
between samples that are combined into one cluster;
-s, --minsamples - (optional parameter, default 2) The number of samples 
(or total weight) in a neighborhood for a point to be considered as a core point. 
This includes the point itself;
-d, --dm - (optional, defaults to None) distance matrix.

-----------------------------------------------------------------------------------
Script mutate.py
-----------------------------------------------------------------------------------

This script is used to convert a consensus sequence in html format to a string 
sequence in which all unreliable positions (with frequent mutations) are marked 
with a special symbol *.

	python3 mutate.py –i consensus.html – o output.fasta –ml “c90 c80”

-i, --input (consensus.html) - html file containing consensus;
-o, --output (output.fasta) - output file name, fasta format;
-ml - levels of nucleotide occurrence, below which nucleotides are noted as mutations.
