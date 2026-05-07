# Quality Control
fastp -i *_R1.fastq.gz -I *_R2.fastq.gz -o *_R1_clean.fastq.gz -O *_R2_clean.fastq.gz -c -h fastp.html -j fastp.json 

# Calculate expression
salmon index -t  transcript.fasta -i SalmonIndex -k 31 --keepDuplicates
salmon quant -i SalmonIndex -l A -1 *_R1_clean.fastq.gz -2 *_R2_clean.fastq.gz -o  *.quant  --validateMappings
