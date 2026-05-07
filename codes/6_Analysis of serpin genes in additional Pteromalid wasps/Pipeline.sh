# Quality Control
fastp -i *_R1.fastq.gz -I *_R2.fastq.gz -o *_R1_clean.fastq.gz -O *_R2_clean.fastq.gz -c -h fastp.html -j fastp.json 

# Align transcriptomic data to the genome
hisat2-build genome.fasta index
hisat2 -p 8 --dta -x index -1 *_R1_clean.fastq.gz -2 *_R2_clean.fastq.gz -S *.sam
samtools sort -@ 8 -o *.bam *.sam 

# Transcript assembly
stringtie -p 8 -o *.gtf *.bam
stringtie --merge -p 8 -o merged.gtf mergelist.txt
gffread merged.gtf -g genome.fasta -w transcript.fa

# Open reading frame prediction was performed according to the procedures described in "1_Processing and annotation of PacBio Iso-Seq data".
# Identification of the serpin3 gene
blastp -query Ppserpin3.ref.fasta -db merged.gtf.pep.fa -evalue 1e-5 -out blastp.out -num_threads 112 -max_target_seqs 1 -outfmt '6 std qlen slen qcovs' 
tblastn -query Ppserpin3.ref.fasta -db genome.fasta -evalue 1e-5 -out tblastn.out -num_threads 112 -outfmt '6 std qlen slen qcovs'

# Calculate expression levels
salmon index -t transcript.fasta -i SalmonIndex -k 31 --keepDuplicates
salmon quant -i SalmonIndex -l A -1 *_R1_clean.fastq.gz -2 *_R2_clean.fastq.gz -o  *.quant --validateMappings
