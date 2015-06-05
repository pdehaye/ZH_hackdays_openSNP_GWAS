"""
    For info, call Paul 
    zero seven six four zero seven five seven nine six
    
    or email 
    paulolivier at gmail
    
    Assumes the following directory structure:
    convert.py
    cleaned_sample/a_23andme_file.txt
    
    Creates the following files:
    selected_genome_XX
    dumped_genome_XX
    parsable_genomes
    snp_count
"""

# Usage:
# python convert.py snps, outputs to snp_count
# counts out snps accross 23andme files

# python convert.py select 80
# selects 23andme files matching all the snps showing up in more than 80% of the files
# generates two files: selected_genome_80 and dumped_genome_80

# python convert.py ped selected_genome_80 snp_count A
# generates two ped files based on slected genomes and snp_count file, adds a column A for the phenotype
# use this twice, concat for the GWAS input to plink


import os, sys

snps_dict = {}

def available_files(directory = "cleaned_sample"):
    """
        Goes through relevant files, according to size and filename
    """
    import os
    return [ filename for filename in list(os.listdir(directory)) if "23andme.txt" in filename\
            and os.stat("cleaned_sample/" + filename).st_size > 22000000]
    
import os

def genome(filename):
    """
        Parses the 23andme files and returns a dict
    """
    with open("cleaned_sample/" + filename, "r") as f:
        genome = {}
        for line in f.xreadlines():
            line = line.strip("\n").strip("\r")
            if line[0] == "#":
                #if "build" in line:
                #    tmp_line = line.split(" ")
                #    print filename, tmp_line[tmp_line.index("build") + 1]
                continue
            #print line[:-1]
            try:
                if "\t" in line:
                    rsid, chromosome, position, genotype = tuple(line.split("\t"))
                else:
                    rsid, chromosome, position, genotype = tuple(line.split(" "))
            except:
                print "     Filename: ", filename
                print "     ", line.split("\t")
            genome[rsid] = genotype
        return genome

def genomes():
    """
        Iterates over genomes in the folder, and yields them. Writes out the parsable ones to a file
    """
    with open("parsable_genomes", "w") as f:
        for filename in available_files():
            try:
                yield filename, genome(filename)
                print "Done", filename
            except:
                print "Skipping " + filename
        f.write(filename+" "+filename+" 1 "+"\n")
    

def select_snps(cutoff_snps = 80):
    """ 
        Selects the SNPs that are present in cutoff percent of the genomes
    """
    cutoff = cutoff_snps/100.
    from collections import Counter
    snp_count = Counter()
    genomes_count = 0
    for filename, genome in genomes():
        for snp in genome:
            snp_count[snp] += 1
        genomes_count += 1
    selected = [(x,y) for x, y in snp_count.items() if y > cutoff*genomes_count]
    dumped = [(x,y) for x, y in snp_count.items() if y <= cutoff*genomes_count]
    #counts.sort(key = lambda (x,y): y, reverse = True)
    print "For indication: "
    print "Parsed ", genomes_count, " genomes"
    print "Dumped ", len(dumped), "snps"
    print "Selected ", len(selected), "snps"
    print "Total", len(dumped) + len(selected), "snps"
    #with open("good_snps_"+str(cutoff_snps), "w") as f:
    #    for snp, count in selected:
    #        f.write(snp+"\n")
    with open("snp_count", "w") as f:
        for x,y in snp_count.items():
            f.write(str(y).rjust(5, "0")+" "+ str(x)+ "\n")
    return [x for (x,y) in selected]

def prune(cutoff_genome = 80):
    """
        Selects genomes that have all the SNPs
    """
    cutoff = cutoff_genome/100.
    counts = {}
    with open("snp_count", "r") as f:
        for line in f.xreadlines():
            count, snp = line.strip("\n").strip("\r").split(" ")
            count = int(count)
            counts[snp] = count
    max_count = max([y for (x,y) in counts.items()])
    good_snps = [x for (x,y) in counts.items() if y > max_count*cutoff]
    with open("selected_genome_"+str(cutoff_genome), "w") as selected_genome_f:
        with open("dumped_genome_"+str(cutoff_genome), "w") as dumped_genome_f:
            for filename, genome in genomes():
                if all(good_snp in genome for good_snp in good_snps):
                    print "Selected ", filename
                    selected_genome_f.write(filename+"\n")
                else:
                    print "Dumped ", filename
                    dumped_genome_f.write(filename+"\n")

def ped():
    pass




if __name__ == "__main__":
    args = sys.argv
    if args[1] == "snps":
        select_snps()
    if args[1] == "select":
        prune(cutoff_genome = int(args[2]))
    if args[1] == "ped":
        ped()
    
    
    