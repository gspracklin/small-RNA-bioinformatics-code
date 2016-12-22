"""
This Script searches for 5' barcodes allowing 1 mismatch. This script assumes
5' barcodes are anchored (starting at index position 0). Once barcode is found
the barcode is removed, as is the last nucleotide (Note: remove this feature if
no complications removing the 3' adaptor).

Name: barcode_splitter.py
Author: George Spracklin
Version: 1.0
"""
import Levenshtein #may need to install, not in standard python package

#Get inputs (fasta, fastq)
usr = input("Enter the file name:")
fasta = open(usr, "r")
z = open("CGTC_barcode", "w")
f = open("AGCG_barcode","w")

#Pre-set barcodes based on Sam Gu 5' barcodes
total = 0
reads_in_file = 0
PP333 = 'AGCG'
PP334 = 'CGTC'

#Search first 4 nucleotides for 5' barcodes, if found write to file
for line in fasta:
    total += 1
    if Levenshtein.distance(line[0:4],PP333) == 1:
        reads_in_file += 1
        f.write(last_line)
        f.write(line[4:-2] + '\n')
    if Levenshtein.distance(line[0:4],PP334) == 1:
        reads_in_file += 1
        z.write(last_line)
        z.write(line[4:-2] + '\n')
    last_line = line

f.close()
z.close()
total = total/2
percentage = (reads_in_file/total)*100
print("the total reads: ", total)
print("reads written to file: ", reads_in_file)
print("percentage of reads contain a barcode:", percentage)
