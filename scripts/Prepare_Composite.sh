mkdir Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome
cd Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome

rename.sh in=Analysis_Dir/Composite_Genome/Ray_Dir/Contigs.fasta out=Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome/contigs.fa prefix=SISRS addprefix=t trd=t

bowtie2-build contigs.fa contigs -p 20
bbmap.sh ref=contigs.fa
samtools faidx contigs.fa

python Analysis_Dir/scripts/Genome_SiteLengths.py Analysis_Dir/Composite_Genome/Ray_Dir/Composite_Genome

