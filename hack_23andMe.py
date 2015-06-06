#!/usr/bin/python

"""
    For info, call Paul 
    zero seven six four zero seven five seven nine six
    
    or email 
    paulolivier at gmail
    
    Assumes the following directory structure:
    convert.py
    a folder with 23andme files
    
    Creates the following files:
    a count file, counting SNPs
    a PED file
    a MAP file
    
"""

# Usage examples:
# 0) help: 
#        ./hack_23andMe.py -h
# 1) generate SNP counts:
#        ./hack_23andMe.py SNP my_23_dir/*
# 2) generates PED file, setting phenotype to 1 for everyone in this file (2 is default), saves to final.{ped,map}:
#        ./hack_23andMe.py --phenotype 1 PED my_23_dir/*
# (one can do twice 2), with different phenotype, and concatenate)
# Run plink:
# ./plink --file final


import os, sys, glob, argparse

def is_valid(filename):
    """Return the opensnp files that are large and have 23andme in their name
    """
    fifteen_mb = 15 * 1000 * 1000
    if os.path.getsize(filename) > fifteen_mb and "23andme" in filename:
        return True
    else:
        print "Warning ", filename, " not large enough"
        return False

def genome(filename):
    """
        Parses the 23andme files and returns a list of (rsid, genotype)
    """
    def process_line(line):
        line = line.strip("\n").strip("\r") 
        try:
            if "\t" in line:
                rsid, chromosome, position, genotype = tuple(line.split("\t"))
            else:
                rsid, chromosome, position, genotype = tuple(line.split(" "))
        except:
            raise Exception(" problematic line: "+ line)
        return rsid, genotype
        
    with open(filename, "r") as f:
        try:
            genome = [process_line(line) for line in f.xreadlines() if line[0] != "#"]
        except Exception as e:
            raise Exception("File "+filename+" has " + str(e))
    return genome

def genomes_iterator(filenames):
    """
        Iterates over the files in filenames, and yields associated genomes
    """
    for filename in filenames:
        if is_valid(filename):
            try:
                yield filename, genome(filename)
                print "Done ", filename
            except Exception as e:
                print "Skipping " + filename + str(e)

def snps(filenames, output ):
    """ 
        Selects the SNPs that are present in cutoff percent of the genomes
    """
    from collections import Counter
    snp_count = Counter()
    genomes_count = 0
    for filename, genome in genomes_iterator(filenames):
        for snp, value in genome:
            snp_count[snp] += 1
        genomes_count += 1
    #counts.sort(key = lambda (x,y): y, reverse = True)
    counts = snp_count.values()
    print "Parsed ", genomes_count, " genomes"
    print "For indication: "
    counts_counter = Counter()
    for count in counts:
        counts_counter[count] += 1
    for repetition, snp_counts in sorted(list(counts_counter.items()), reverse = True):
        print "        ", snp_counts, " times a snp appearing ", repetition, " times"
    print "Writing ", output
    with open(output, "w") as f:
        for x,y in sorted(list(snp_count.items()), reverse = True, key = lambda (x,y) : y):
            f.write(str(y).rjust(5, "0")+" "+ str(x)+ "\n")


def good_snps_reader(snp_count_file, cutoff_genome, output):
    """
        Reads in good SNPs and creates MAP file
    """
    cutoff = cutoff_genome/100.
    with open(snp_count_file, "r") as f:
        counts = map(lambda (count, snp) : (snp, int(count)),  \
            [line.strip("\n").strip("\r").split(" ") for line in f.xreadlines()])
    max_count = max([y for (x,y) in counts])
    good_snps = sorted([x for (x,y) in counts if y > max_count*cutoff])
    print "Selected at ", cutoff_genome, " percent ", len(good_snps), " snps out of ", len(counts)
    
    map_output = output.split(".")[0]+".map"
    with open(map_output, "w") as f:
        for snp in good_snps:
            f.write("0 "+snp+ " 0 1000\n")
    
    return good_snps

def convert_allele(allele):
    """
            AA ---> A A
            -A ---> 0 A
            -- ---> 0 0
    """
    tmp = allele.replace("-", "0")
    if len(tmp) == 1:
        tmp = 2 * tmp 
    tmp = " ".join(list(tmp))
    return tmp

def ped(snp_count_file, cutoff_genome, phenotype_value, genome_filenames, output):
    """
        Selects genomes that have all the SNPs, creates standard PED file 
        
    """
    
    
    phenotype_value = str(phenotype_value)
    assert phenotype_value in ["2", "1"] # Expects 2, 1 for affected/unaffected
    good_snps = good_snps_reader(snp_count_file, cutoff_genome, output)

    selected, rejected = 0, 0
    with open(output, "w") as f:
            for filename, genome in genomes_iterator(genome_filenames):
                genome_dict = dict(genome)
                id = filename.split(os.path.sep)[-1].strip("23andme.txt")
                try:
                    print "Selected ", filename
                    v = [convert_allele(genome_dict[good_snp]) for good_snp in good_snps]
                    selected += 1
                    f.write(" ".join([id, id, "0", "0", "0", phenotype_value])+ "    " + "  ".join(v) + "\n")
                    # family ID, individual ID, Paternal ID, Maternal ID, sex, phenotype
                except:
                    print "Rejected", filename
                    # This guy is missing a good snp
                    rejected += 1

    print "Rejected: ", rejected
    print "Selected: ", selected

    # I decided to create complex PED file so command line is simple
    # THESE WOULD BE MORE COMPLEX OPTIONS TO USE: --compound-genotypes --no-parents --no-sex --no-fid --missing-genotype - 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Manipules 23andMe files.')
    parser.add_argument('action', help = "Main action to perform. SNP creates an SNP count file, PED creates PED and MAP files", choices= ["SNP", "PED"])
    parser.add_argument('--cutoff', type = int, nargs = 1, default = 90, help = 'Only consider SNPs that are present in this percent of genomes, only necessary for PED action')
    parser.add_argument('--phenotype', default = 2, help = "the phenotype value to write to the final output; 2 is disease present, 1 is absent", choices= ["1", "2"])
    parser.add_argument('--snpoutput', nargs = 1, default = 'SNPfile.txt', help = "The file to use as for output for SNP count procedure")
    parser.add_argument('genomes', nargs = "*", help = "The 23andMe files to consider")
    parser.add_argument('--pedoutput', default = 'output.ped', help = "The base filename to use as for output for PED procedure, for .ped and .map files")
    

    args = parser.parse_args()
    if args.action == "SNP":
        snps(args.genomes, args.snpoutput)
    if args.action == "PED":
        ped(args.snpoutput, args.cutoff, args.phenotype, args.genomes, args.pedoutput)
    
    
    