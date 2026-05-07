# Quality Control
fastp -i *_R1.fastq.gz -I *_R2.fastq.gz -o *_R1_clean.fastq.gz -O *_R2_clean.fastq.gz -c -h fastp.html -j fastp.json 

# Calculate expression
salmon index -t  transcript.fasta -i SalmonIndex -k 31 --keepDuplicates
salmon quant -i SalmonIndex -l A -1 *_R1_clean.fastq.gz -2 *_R2_clean.fastq.gz -o  *.quant  --validateMappings

# Obtain Bigwig files for visualization in IGV
# The bigWig files were generated following the workflow provided in PipeRNAseq
# merge bigwig files
bigWigMerge PpCAa.Ppup.sorted.bedGraph.bw PpCAb.Ppup.sorted.bedGraph.bw PpCAc.Ppup.sorted.bedGraph.bw PpCA.merge.bedGraph
sort -k1,1 -k2,2n PpCA.merge.bedGraph > PpCA.merge.sort.bedGraph
bedGraphToBigWig PpCA.merge.sort.bedGraph Ppup.ChromInfo.txt PpCA.merge.bw
