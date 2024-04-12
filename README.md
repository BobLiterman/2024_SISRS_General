# Running SISRS for Diagnostic SNP Development
#### Dr. Robert Literman
#### Center for Food Safety and Applied Nutrition
#### US Food and Drug Administration
---
This repo is designed to provide a more manual walkthrough for using SISRS to identify diagnostic SNPs. If you clone the repo, all the scripts and preferred folder structure will be automatically in place. 

### Step 1: Read Data Organization
- **Software Requirements**: Your favorite read trimmer
- I prefer to use this helpful naming convention for trimmed reads: For samples of *Homo sapiens*, I name the trimmed reads HomSap_SIMPLE_ID. This way, if you want to grab reads from a specific species you can do simple soft-links (See below)

- It is not required, but for ease trimmed reads from each species should go into separate directories

### Step 2: Composite Genome Assembly

- **Software Requirements**: Ray, mpirun, python, samtools v.1.3.1+, bowtie2, bbtools

1. If your reads are **NOT** organized by species, and are sensibly named, you can leave your data where it is and just use softlinks as such:
   
        ```
        # Make and enter directory
        mkdir Analysis_Dir/Reads/Species_Reads
        cd Analysis_Dir/Reads/Species_Reads

        # Create a folder for each species
        mkdir HomSap GorGor PanPan

        # Softlink reads into each folder
        for f in *; do cd $f; cp -as 'Analysis_Dir/Reads/Trimmed_Reads/'$f*'gz' .;cd ..;done
        ```

2. To generate a composite read dataset where each species is represented equally, use [Read_Subsetter.py](scripts/Read_Subsetter.py) (10X coverage has always worked well in my hands):
    ```
    # python Read_Subsetter.py -g GROUP_GENOME_SIZE_ESTIMATE -c COVERAGE -r SPECIES_READ_DIR -o OUTPUT_DIR

    # To generate 10X coverage for a 4.5Gb genome:
    python Read_Subsetter.py -g 4500000000 -c 10 -r Analysis_Dir/Reads/Species_Reads -o Analysis_Dir/Reads/Composite_Reads

    # To generate 20X coverage for a 1Mb genome:
    python Read_Subsetter.py -g 1000000 -c 20 -r Analysis_Dir/Reads/Species_Reads -o Analysis_Dir/Reads/Composite_Reads

    ```
   - **NOTE:** The read subsetting step always outputs .fq.gz files with a _1/_2 naming convention. It assumes input files are .fastq.gz, where forward and reverse reads are _1.fastq.gz/_2.fastq.gz. **If your reads are not _1.fastq.gz/_2.fastq.gz**, convert them or make the relevant changes in Read_Subsetter.py. 
    
        ```
        # If your input reads are _R1.fastq/_R2.fastq, you could run:

        sed -i 's/_1.fastq.gz/_R1.fastq/g' Read_Subsetter.py
        sed -i 's/_2.fastq.gz/_R2.fastq/g' Read_Subsetter.py
        sed -i 's/fastq.gz/fastq/g' Read_Subsetter.py
        ```

3. Assemble Composite Genome
    ```
    # Load modules if needed

    # module load openmpi
    # module load Ray

    # E.g., if you have 10 nodes with 20 CPU per node, you could run:

    mpirun -np 200 Ray -detect-sequence-files Analysis_Dir/Reads/Composite_Reads -o Analysis_Dir/Composite_Genome/Ray_SOMEUSEFULID
    ```
    - **Note**: Ray output folders must start with 'Ray_'

4. Prepare composite genome
   ```
    # Load modules if needed

    # module load bbmap
    # module load samtools
    # module load bowtie2
   
    mkdir Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome
    cd Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome

    rename.sh in=Analysis_Dir/Composite_Genome/Ray_Dir/Contigs.fasta out=Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome/contigs.fa prefix=SISRS addprefix=t trd=t

    bowtie2-build contigs.fa contigs -p PROCESSOR_COUNT
    bbmap.sh ref=contigs.fa
    samtools faidx contigs.fa

    python Analysis_Dir/scripts/Genome_SiteLengths.py Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome
   ```
