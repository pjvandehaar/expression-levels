This program implements a simplistic approach to checking expression levels of genes.

It takes a Fasta file of RNAseq reads and a csv file of genes. It compares every read with every gene and computes a match score.

The match scores  computed using local alignment with an infinite gap penalty. So, they are the complements of hamming distance, aligned between the reads.

Unexpressed genes will just show a biased Binomial distribution. Expressed genes will have a few really good matches.

## Usage
1. Get some Fasta files with the extension `.fas` in the working directory.
2. Put some gene names and sequences into the `query_sequences.csv` file.
3. Run `./get_match_scores.py` under Python2 or Python3.
4. Sit and watch the output for a few years.
