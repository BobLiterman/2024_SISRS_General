#!/usr/bin/env python3

# TRIMMED READS SHOULD BE QC CHECKED PRIOR TO THIS STEP TO REMOVE DATASETS WITH LOW QUALITY OR HIGH DUPLICATION RATES
# This script performs read subsetting for a SISRS composite genome assembly
# This script calls reformat.sh (part of BBMap Suite), which must be installed and in your path
# All reads for all taxa should be in .fastq.gz format (To change this, find/replace this script, replacing '.fastq.gz' with your chosen extension)
# Paired-end read files must be identically basenamed and end in _1.fastq.gz/_2.fastq.gz
#
# Arguments: 
# (1) -g/--genomesize (REQUIRED): An estimate of the  group genome size (e.g. 3500000000 for a primate dataset)
# (1) -c/--coverage (OPTIONAL): When sampling reads, how many X coverage total do you want? (e.g. 10X for primates = 35Gb total)
# (1) -r/--reads (OPTIONAL): Directory containing folders of reads [Default: ../Reads/TrimReads]
# (1) -o/--output (OPTIONAL): Directory to output subset reads [Default: ../Reads/SubsetReads]
#
# Output: For an analysis with x taxa groups, a genome size estimate of g, and a final desired coverage of c, this script will subset reads from each taxon group down to (c*g)/x bases

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call
import pandas as pd
import re
import argparse
from pathlib import Path

#Set cwd to script location
script_dir = os.path.dirname(os.path.abspath(__file__))

#Set Genome size estimate, desired coverage, data location
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-g','--genomesize',action='store',required=True,type=int)
my_parser.add_argument('-c','--coverage',action='store',default=10,type=float)
my_parser.add_argument('-r','--reads',action='store',default=path.dirname(path.abspath(script_dir))+"/Reads/TrimReads",type=str)
my_parser.add_argument('-o','--output',action='store',default=path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads",type=str)
args = my_parser.parse_args()

genomeSize = int(args.genomesize)
genomeCov = float(args.coverage)
readPath = os.path.normpath(str(args.reads))
outPath = os.path.normpath(str(args.output))

#Check/set path where trimmed read files live
if not path.isdir(readPath):
    sys.exit("Path to reads does not exist.")
trim_read_dir = readPath

#Check if output folder exists and is empty
if path.isdir(outPath):
    if len(os.listdir(outPath)) > 0:
        sys.exit("Output path not empty.")
else:
    os.mkdir(outPath)
subset_output_dir = outPath

#Create folder for subset output logs
os.mkdir(outPath+"/subsetOutput")
subset_log_dir = outPath+"/subsetOutput/"

#Set taxa folders to only those containing .fastq.gz
trim_read_tax_dirs = sorted(list(set([os.path.dirname(f) for f in list(Path(trim_read_dir).rglob("*.fastq.gz" ))])))

if len(trim_read_tax_dirs) == 0:
    sys.exit("No subfolders of read directory contain .fastq.gz files.")

#Calculate subset depth
subsetDepth = int((genomeCov*genomeSize)/len(trim_read_tax_dirs))

print("Based on a genome size estimate of " + str(genomeSize) + " bp, and with " + str(len(trim_read_tax_dirs)) + " species/groups, the requested subset depth at " + str(genomeCov) + "X coverage is " + str(subsetDepth) + " bp per species/group")

#Initialize Pandas DF to get base counts and lists of paired and single-end files
df = pd.DataFrame(columns=['Taxon','Dataset','Basecount'],dtype=object)
compiled_paired = list()
compiled_single_end = list()

#For each taxa directory...
for tax_dir in trim_read_tax_dirs:
    #List all files and set output dir
    taxon_ID = path.basename(tax_dir)

    left_pairs = list()
    right_pairs = list()
    single_end = list()
    taxon_list = list()
    dataset_list = list()
    basecount_list = list()
    
    #Find fastq.gz files 
    files = sorted(glob(tax_dir+"/*.fastq.gz"))

    #Find fastq.gz files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_1.fastq.gz" in s]
    right_files = [s for s in files if "_2.fastq.gz" in s]

    #Strip _1_Trim.fastq.gz/_2_Trim.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_2.fastq.gz', '') for x in right_files]
    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    for pair in paired_files:
        left_pairs.append(pair+"_1.fastq.gz")
        right_pairs.append(pair+"_2.fastq.gz")
    paired_files = sorted(left_pairs + right_pairs)

    single_end = [x for x in files if x not in paired_files]

    compiled_paired = compiled_paired + paired_files
    compiled_single_end = compiled_single_end + single_end
    
    #Remove .fastq.gz from lists to make naming easier
    left_pairs = [x.replace('_1.fastq.gz', '') for x in left_pairs]
    right_pairs = [x.replace('_2.fastq.gz', '') for x in right_pairs]
    single_end = [x.replace('.fastq.gz', '') for x in single_end]

    #Count bases in single-end files if present...
    if len(single_end) > 0:
        for x in single_end:
            count_command = [
                'reformat.sh',
                'in={}'.format(x+'.fastq.gz')]
            proc = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            o, e = proc.communicate()

            #Munge output to get base count
            count = re.search('HERE (.*) HERE_2', str(e).replace('(100.00%) \\t','HERE ').replace(' bases (100.00%)\\n\\nTime:',' HERE_2'))
            baseCount = count.group(1)
            taxon_list.append(taxon_ID)
            dataset_list.append(path.basename(x))
            basecount_list.append(baseCount)

    #Count bases in paired-end files if present...
    if(len(left_pairs) == len(right_pairs) & len(left_pairs) > 0):
        for x in range(len(left_pairs)):
            dataset_list.append(path.basename(left_pairs[x]))
            count_command = [
                'reformat.sh',
                'in={}'.format(left_pairs[x]+'_1.fastq.gz'),
                'in2={}'.format(right_pairs[x]+'_2.fastq.gz')]
            proc = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            o, e = proc.communicate()

            count = re.search('HERE (.*) HERE_2', str(e).replace('(100.00%) \\t','HERE ').replace(' bases (100.00%)\\n\\nTime:',' HERE_2'))
            baseCount = count.group(1)
            basecount_list.append(int(baseCount))
            taxon_list.append(taxon_ID)
    list_of_tuples = list(zip(taxon_list,dataset_list, basecount_list))
    df = df.append((pd.DataFrame(list_of_tuples, columns = ['Taxon','Dataset', 'Basecount'])))

df = df.reset_index(drop=True)
df["Basecount"] = pd.to_numeric(df["Basecount"])
df = df.sort_values(by=['Taxon'])
df["ToSample"] = 0
taxa_list = sorted(df.Taxon.unique())

for taxa in taxa_list:
    taxaDF = df[df['Taxon']==taxa]

    #Check if taxa has enough total coverage. If not, use all reads and warn user.
    if(taxaDF["Basecount"].sum() < subsetDepth):
        print("NOTE: " + taxa + " samples only have a total of " + str(taxaDF["Basecount"].sum()) + " bp. All " + taxa + " reads will be used.\nWARNING: Having under the required read threshold will have negative consequences on composite genome assembly and site calling.\n")
        df.loc[df["Taxon"] == taxa, ["ToSample"]] = list(df[df["Taxon"] == taxa].Basecount)

    elif(taxaDF["Basecount"].sum() >= subsetDepth):
        #If there is enough total data, check sample values to optimize sampling
        dataset_count = taxaDF.shape[0] #Count samples for taxa
        per_dataset = subsetDepth/dataset_count #Calculate bases needed per sample if evenly distributed
        s = pd.Series(taxaDF["Basecount"]) #Grab basecount colulmn

        if((s < per_dataset).any()):
            subset_countdown =  subsetDepth #Create countdown variable

            while((s < per_dataset).any()):

                low_datasets = list(taxaDF[taxaDF.Basecount<per_dataset]["Dataset"]) #Find datasets with too few bases
                high_datasets = list(taxaDF[taxaDF.Basecount>=per_dataset]["Dataset"]) #Find datasets with enough bases (currently)

                sum_low = taxaDF[taxaDF.Basecount<per_dataset].Basecount.sum() #Sum the bases of the low datasets

                print("Not all " + taxa + " samples have enough for even read sampling (" + ",".join(low_datasets) + "). Will take all reads from data poor samples and compensate from larger samples...")
                df.loc[df.Dataset.isin(low_datasets), ["ToSample"]] = list(df[df.Dataset.isin(low_datasets)].Basecount) #Set low datasets to use all reads

                taxaDF = taxaDF[taxaDF.Dataset.isin(high_datasets)] #Reset taxa DF to remove low samples
                dataset_count = taxaDF.shape[0] #Count remaining datasets
                subset_countdown = subset_countdown - sum_low

                per_dataset = int(subset_countdown/dataset_count) #Reset per_dataset

                s = pd.Series(taxaDF.Basecount) #Recount basecounts for while loop

            df.loc[df.Dataset.isin(high_datasets), ["ToSample"]] = per_dataset
        else:
            df.loc[df.Taxon == taxa, ["ToSample"]] = per_dataset

df.to_csv(subset_output_dir+"/Subset_Scheme.csv",index=False)

compiled_paired = [path.basename(x).replace('_1.fastq.gz', '') for x in compiled_paired]
compiled_paired = [path.basename(x).replace('_2.fastq.gz', '') for x in compiled_paired]
compiled_single_end = [path.basename(x).replace('.fastq.gz', '') for x in compiled_single_end]

#Subset bases in single-end files if present...
if len(compiled_single_end) > 0:
    single_end_DF = df[df.Dataset.isin(compiled_single_end)]
    for row in single_end_DF.itertuples():
        subset_command = [
            'reformat.sh',
            'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+".fastq.gz"),
            'out={}'.format(subset_output_dir+"/"+str(row.Dataset)+'_GenomeReads.fastq.gz'),
            'samplebasestarget={}'.format(int(row.ToSample)),
            'ow=t',
            '&>',
            '{outDir}out_{fileName}_Subset'.format(outDir=subset_log_dir,fileName=row.Dataset)]
        check_call(subset_command)

#Subset bases in paired-end files if present...
if len(compiled_paired) > 0:
    paired_DF = df[df.Dataset.isin(compiled_paired)]
    for row in paired_DF.itertuples():
        subset_command = [
            'reformat.sh',
            'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_1.fastq.gz"),
            'in2={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_2.fastq.gz"),
            'out={}'.format(subset_output_dir+"/"+str(row.Dataset)+'_GenomeReads_1.fq.gz'),
            'out2={}'.format(subset_output_dir+"/"+str(row.Dataset)+'_GenomeReads_2.fq.gz'),
            'samplebasestarget={}'.format(int(row.ToSample)),
            'ow=t',
            '&>',
            '{outDir}out_{fileName}_Subset'.format(outDir=subset_log_dir,fileName=row.Dataset)]
        check_call(subset_command)

