"""
This script takes SAM file and parses based on chromosme and base position.
Creates two files, sense reads and antisense reads. Then plots data using
matlibplot. (Note: Matlibplot settings may need changing). 

Script:nucleotide_density_plot.py
Author: George Spracklin
v1.0
"""

#inputs
usr0 = input("What file (SAM format only) :")
f = open(usr0, 'r')
usr3 = input("What chromosome? :")
usr1 = int(input("What is the starting position? :"))
usr2 = int(input("What is the ending position? :"))
rpm = float(input("What is the total number of reads (in Millions)? :"))
norm_factor = float(input("What is the normalization fraction? :"))
import matplotlib.pyplot as plt
import matplotlib
from collections import Counter
import numpy as np

##############################################
#find reads in SAM file using Chr # and position number
##############################################
mapped_siRNAs = []
total_reads = 0
for line in (f):
    parts = line.split()
    if len(parts) > 4:
        if parts[2] == usr3: #Chromosome
            if int(parts[3]) >= usr1 and int(parts[3]) < usr2: #basepair number
                mapped_siRNAs.append(parts)
print("finished scanning SAM file")

##############################################
#Separates Sense and Antisense reads
##############################################
antisense_reads = []
sense_reads = []
anti_count = 0
for line in mapped_siRNAs:
    read = len(line[9])
    start = int(line[3])
    if line[1] == '16':
        antisense_reads.append((start,read))
        anti_count += 1
    if line[1] == '0':
        sense_reads.append((start,read))
print("finished sorting sense/antisense")
print("total mapped reads equals: ", anti_count)

##############################################
#Sorting Sense and Antisense Reads
##############################################
#Antisense Reads
positions_antisense = []
for length in antisense_reads:
    counter = length[0]
    for item in range(length[1]):
        positions_antisense.append(counter)
        counter +=1
positions_antisense.sort()
A = Counter(positions_antisense)
#s = open("antisense_reads", 'w')
#for k,v in  A.most_common():
#    s.write( "{} {}\n".format(k,v) )

#Sense Reads
positions_sense = []
for length in sense_reads:
    counter = length[0]
    for item in range(length[1]):
        positions_sense.append(counter)
        counter +=1
positions_sense.sort()
S = Counter(positions_sense)
#q = open("sense_reads", 'w')
#for k,v in  S.most_common():
#    q.write( "{} {}\n".format(k,v) )

##############################################
#Normalization
##############################################
f_count = 0
for key, value in A.items():
    f_count += 1
    A[key] = ((value / rpm)/ norm_factor)

for key, value in S.items():
    S[key] = ((value / rpm)/ norm_factor)

print("finished normalizing RPM/norm_factor")

##############################################
#Create histogram
##############################################
N = usr2 - usr1

#Create a Control_dict
Control_dict = {}
for k in range(N):
    k_new = k + usr1
    Control_dict[k_new] = 0


#Create key for any missing nucleotides and set value 0
for k in Control_dict:
    A.setdefault(k, 0)
    S.setdefault(k, 0)

#Transform dictionary to a list
antisense_x_data = []
antisense_y_data = []
for k, v in A.items():
    antisense_x_data.append(k)
    antisense_y_data.append(v)

sense_x_data = []
sense_y_data = []
for k, v in S.items():
    sense_x_data.append(k)
    sense_y_data.append(v)
print("finished merging dictionaries")
#plot data
#fig = plt.figure(1)
f, (ax1, ax2) = plt.subplots(2, sharex='col', sharey='row')

ax1.plot(antisense_x_data, antisense_y_data, color ='blue')
ax1.set_title('Antisense', fontsize=20)
ax1.set_yscale("log", basey=10)
ax1.set_ylim(0.5,3000)

#ax2 = plt.subplot(212)
ax2.plot(sense_x_data, sense_y_data, color ='red')
ax2.set_title('Sense', fontsize=20)
ax2.set_yscale("log", basey=10)
ax2.set_ylim(0.5,3000)

# Set common labels
f.text(0.5, 0.04, 'genomic locus', ha='center', va='center')
f.text(0.06, 0.5, 'siRNA coverage (log scale)', ha='center', va='center', rotation='vertical')
f.savefig("%s.log.pdf" % usr0)
