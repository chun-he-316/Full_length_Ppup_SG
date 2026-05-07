1. IsoSeq3 pipeline
# merge two datasets and run together for the workflow. This will make a consensus of the subreads:
pbmerge -o CaSGMerged.subreads.bam split.5p--bc1002_3p.subreads.bam  split.5p--bc1001_3p.subreads.bam
# Step1, CCS generation
ccs CaSGMerged.subreads.bam CaSGMerged.ccs.bam --min-rq 0.9
# Step 2 Primer removal and demultiplexing
lima CaSGMerged.ccs.bam primer.fasta CaSGMerged.ccs.fl.bam --isoseq --peek-guess
# Step3, refine
isoseq3 refine CaSGMerged.ccs.fl.NEB_5p--NEB_Clontech_3p.bam primer.fasta CaSGMerged.ccs.fl.flnc.bam --require-polya
# Step4, cluster
isoseq3 cluster CaSGMerged.ccs.fl.flnc.bam CaSGMerged.ccs.fl.flnc.finalcluster.bam --verbose --use-qvs

2. ToFU Cupcake collapse
minimap2 -ax splice -t 30 -uf --secondary=no -C5 JXF_scaffold_V1.1.rmBac.fasta CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta > CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sam
sort -k 3,3 -k 4,4n CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sam > CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted.sam
python ./cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py --input CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta -s CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted.sam --dun-merge-5-shorter -o CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted
# This step will generate read_stat.txt and abundance.txt files.
python ./cDNA_Cupcake/cupcake/tofu/get_abundance_post_collapse.py CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted.collapsed CaSGMerged.ccs.fl.flnc.finalcluster.cluster_report.csv
# Then filter 5end degradation reads
python ./cDNA_Cupcake/cupcake/tofu/filter_away_subset.py CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted.collapsed

3. SQANTI3
python ./SQANTI3/sqanti3_qc.py CaSGMerged.ccs.fl.flnc.finalcluster.hq.fasta.Ppup.sorted.collapsed.filtered.gff  Ppup.RefSeq.reduced.bed12.geneid.gtf JXF_scaffold_V1.1.rmBac.fasta  -d ./ --fl_count CaSGMerged.cupcake.pbIsoCollapse.abundance.txt -t 5 -n 10 -o sqanti3 
python ./SQANTI3/sqanti3_filter.py rules sqanti3_classification.txt --gtf sqanti3_corrected.gtf -o sqanti3_filter -d ./ 

python ./software/gtfmerge/gtfmerge.py  --output_prefix  Ppup.final.gtf  --tol_tts 500 --tol_tss 500 Ppup.OGS2.1.gtf  sqanti3_filter.filtered.addx.gtf

4. TransDecoder
perl ./TransDecoder-master/util/gtf_genome_to_cdna_fasta.pl Ppup.final.gtf JXF_scaffold_V1.1.rmBac.fasta > transcripts.fasta
perl ./TransDecoder-master/util/gtf_to_alignment_gff3.pl Ppup.final.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta 
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /data/public/database/uniprot/uniprot  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 112 > blastp.outfmt6 
hmmsearch --cpu 112 --domtblout pfam.domtblout Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep 
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6  --single_best_only
