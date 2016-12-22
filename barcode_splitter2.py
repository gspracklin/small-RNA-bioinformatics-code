"""
This Script searches for 5' barcodes allowing 1 mismatch. This script assumes
5' barcodes are anchored (starting at index position 0). Once barcode is found
the barcode is removed, as is the last nucleotide (Note: remove this feature if
no complications removing the 3' adaptor).

Name: barcode_splitter.py
Author: George Spracklin
Version: 1.1
"""
import Levenshtein #may need to install, not in standard python package

#Get inputs (fasta, fastq)
usr = input("Enter the file name:")
fasta = open(usr, "r")
barcode = input("Enter bacode: ")
f = open(usr + barcode,"w")

reads_in_file = 0
correct = 2
for line in fasta:
    if '@' in line:
        header = line
        correct = 2
    if correct == 1:
        f.write(line)
    if Levenshtein.distance(line[0:4],barcode) == 1: #Search first 4 nucleotides for 5' barcodes, if found write to file
        reads_in_file += 1
        f.write(header)
        f.write(line)
        correct = 1

f.close()
print("the total reads: ", reads_in_file)
