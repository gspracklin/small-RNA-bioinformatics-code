"""
This script writes and counts all reads mapping to user defined chromosome interval. 

v1.0
Author:George Spracklin
"""


usr0 = input("What file (SAM format only) :")
f = open(usr0, 'r')
w = open(usr0 + '_oma-1', 'w')
usr3 = input("What chromosome? :")
usr1 = int(input("What is the starting position? :"))
usr2 = int(input("What is the ending position? :"))


total_reads = 0
for line in f:
    total_reads += 1
    parts = line.split()
    if len(parts) > 4: #removes header from analysis
        if len(parts[9]) < 30: # analyzes reads short than 30bp
            if parts[2] == usr3: #Chromosome selection
                if int(parts[3]) >= usr1 and int(parts[3]) < usr2: #basepair number
                    w.write(line)
print("Total reads,", total_reads)
