import re
import sys
import os.path
import csv


# -------------------------------------- Subroutines ----------------------------------------- #
def help(error_no):
    print "Give new IDs to gene IDs for CORB and LUC projects"
    print
    print "Usage:"
    print "    %s gff3_file FNUM output_file" % os.path.basename(sys.argv[0])
    print
    print "Acronyms:"
    print "    FNUM: First number of the gene desired for each chromosome"
    sys.exit(error_no)

def add2dict(key, value, myDict):
    try:
        if value not in myDict[key]:    # Add value if it does not exist
            myDict[key].append(value)
    except KeyError:
        myDict[key] = [value]


def FindAndReplace(fname, search, replace, oname):        
    search = re.sub(r"\.", r"\\.", search)   # If string has a dot add '\' before to correct finding
    out = open(oname, 'a')
    with open(fname, "r") as dataf:
        reader = csv.reader(dataf, delimiter="\t")
        c = 0                   # Just for control. At least one match must be found
        for line in reader:
            # if search in line[8]:
            if re.search(r'\b' + search + r'\b', line[8]):
                c = 1
                if 'Coming_from' in line[8]:
                    f8_1, f8_2 = line[8].split("Coming_from=", 1)
                    line[8] = re.sub(search, replace, f8_1) + 'Coming_from=' + f8_2
                    # print line[8]
                    # line[8] = f8_1.replace(search, replace) + 'Coming_from=' + f8_2
                else:
                    line[8] = re.sub(search, replace, line[8])
                new_line = '\t'.join(line) + '\n'
                out.write(new_line)
        out.close()
        if c == 0:
            print "Error. No matching IDs."
            sys.exit()


# ----------------------------------------- Arguments ----------------------------------------- #
if len(sys.argv) != 4: help(1)
gff3_file = sys.argv[1]
fnum = int(sys.argv[2])
output_file = sys.argv[3]

if not os.path.exists(gff3_file):
    print >> sys.stderr, "GFF3 file (%s) was not found!" % gff3_file
    sys.exit(1)
if os.path.exists(output_file):
    os.remove(output_file)

# ----------------------------------------- First step ----------------------------------------- #
# Process gff3 input file
# Save all gene IDs into a dictionary
old = {}
# Read input gff3
print "1. Reading input file ..."
with open(gff3_file, "r") as gff3:
    reader = csv.reader(gff3, delimiter="\t")
    for gff3_line in reader:
        if gff3_line[0] == "#":
            continue
        else:
            if gff3_line[2] == "transcript":
                parent = gff3_line[8].split(";", 1)[0]
                gene_id = parent.replace("Parent=", "")
                # chr_num = gene_id.split("g", 1)[0].replace("Solyc", "").replace("CORB", "")
                # chr_num = 'chr_' + chr_num
                add2dict(gff3_line[0], gene_id, old)
            else:
                continue    # go to next line

# ----------------------------------------- Second step ----------------------------------------- #
# Rename gene IDs and create a new gff3 file with the new IDs
print "2. Renaming genes ..."
new = {}
for chrom, genes in old.iteritems():
    # print chrom, "=>", genes
    chr_num = chrom.replace("RC-super_SL2.50ch", "")
    new_num = fnum
    for g in genes:
        if g.endswith('C') or g.endswith('C_bis'):
            new_id = 'SolycLC' + str(chr_num) + 'g' + str(new_num).zfill(6) + 'C'    # str(new_num).zfill(6) --> In order to mantain the number of digits independently
        else:
            new_id = 'SolycLC' + str(chr_num) + 'g' + str(new_num).zfill(6)
        # Find OLD gene ID and replace with NEW one
        FindAndReplace(gff3_file, g, new_id, output_file)
        add2dict(chrom, new_id, new)
        new_num += 10


# ----------------------------------------- Third step ----------------------------------------- #
# Create a conversion file with two columns: old_gene_ID    new_gene_ID
print "3. Creating conversion table ..."
try:
    os.remove('conversion_table.txt')
except OSError:
    pass
conv_file = open('conversion_table.txt', 'a')
for (k1,v1), (k2,v2) in zip(old.items(), new.items()):
    for (s1), (s2) in zip(v1, v2):
        conv_file.write(s1 + "\t" + s2 + "\n")
conv_file.close()


# for key,val in old.items():
#     print key, "=>", val
# for key,val in new.items():
#     print key, "=>", val
