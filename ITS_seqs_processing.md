Workflow followed when processing fungal ITS region sequencing data
- raw sequences preprocess
- OTU construction and taxonomy assignment
- preparing for downstream analyses


## RAW SEQUENCES PREPROCESS


#### DEMULTIPLEXING 

~~~
split_libraries_fastq.py -i raw_its_r1.fastq -o . -m map_gard2.txt -b out.barcodes.fastq --barcode_type 12 --rev_comp_mapping_barcodes --store_demultiplexed_fastq
~~~

#### QUALITY FILTERING
min. qual mean 25 (first trimmed right end at 25 quality), no N bases
~~~
prinseq-lite.pl -fastq gard.fastq -min_qual_mean 25 -ns_max_n 0 -trim_qual_right 25 -rm_header -out_format 4 -out_good gard_q25 -out_bad null -graph_data gard_q25.gd
~~~

#### ITS extraction
~~~
ITSx -i ../gard_q25.fasta -o gard --preserve T --cpu 12
~~~
> Number of sequences in input file:              1742893
Sequences detected as ITS by ITSx:      1630014


#### LENGTH FILTERING & TRIMMING
min length 150
~~~
prinseq-lite.pl -fasta itsx/gard.ITS1.fasta -out_good gard_l150 - out_bad null -min_len 150 -line_width 0
~~~
> Input sequences: 1,294,433
>-	Good sequences: 991,961 (76.63%)

Check filtering success
~~~
prinseq-lite.pl -fasta rum.ITS1_l150.fasta -stats_len -out_good null -out_bad null
~~~
>- stats_len       max     221
>-  stats_len       mean    174.64
>- stats_len       median  172
>- stats_len       min     150
>- stats_len       mode    163
>- stats_len       modeval 145498
>- stats_len       range   72
>- stats_len       stddev  13.34


## OTU CONSTRUCTION
USEARCH needs “.” instead of “_” in the sequence header
~~~
sed 's/_/./g' ../gard_l150.fasta > seqs_dot.fna
~~~
Dereplication
~~~
usearch81_64_new -derep_fulllength seqs_dot.fna -fastaout uniq.fna -sizeout
~~~
>-	933570 seqs, 161613 uniques, 130761 singletons (80.9%)
>-	Min size 1, median 1, max 62581, avg 5.78
	
OTU clustering (98.5% similarity), singletons removed
~~~
usearch81_64_new -cluster_otus uniq.fna -minsize 2 -otus otus.fna -relabel OTU -otu_radius_pct 3
~~~
>	2526 OTUs, 409 chimeras (1.3%)

## TAXONOMY ASSIGNMENT
Taxonomy assignment; Qiime - BLAST against UNITE 7.2
~~~
parallel_assign_taxonomy_blast.py -i ../otus.fna -o. -r /mnt/data/chomic/databases/unite/unite_qiime_17_12_01/sh_refs_qiime_ver7_dynamic_s_01.12.2017.fasta  -O 16 -t /mnt/data/chomic/databases/unite/unite_qiime_17_12_01/sh_taxonomy_qiime_ver7_dynamic_s_01.12.2017.txt
~~~

Convert BLAST output to UTAX compatible
>get rid of extra \n  within sequences in otus.fna
~~~
prinseq-lite.pl -fasta otus.fna -out_bad null -out_good otus_line -line_width 0
~~~
Within BLAST output, sort OTUs ascending
~~~
sed 's/OTU//1' blast/otus_tax_assignments.txt | sort -g | sed 's/^/OTU/' > blast/otus_tax_assignments_sorted.txt
~~~
Add empty rows to merge with otus.fna, get rid of OTU number
~~~
sed -e 'G' blast/otus_tax_assignments_sorted.txt | cut -f 2 > blast/blast_taxonomy.txt
~~~
Reformat taxonomy syntax for utax; delete "Ambigous_taxa" fields
~~~
sed -e 's/;Ambiguous_taxa//g; s/D_0__/tax=d:/g; s/__/:/g; s/;/,/g; s/D_1/p/g; s/D_2/c/g; s/D_3/o/g; s/D_4/f/; s/D_5/g/g; s/D_6/s/g; s/ /_/g; s/tax=d/;tax=d/g' blast/blast_taxonomy.txt > blast/blast_taxonomy_utax.txt

~~~
Merge OTU sequences and UTAX comaptible taxonomy
~~~
paste -d "" otus_line.fasta blast_taxonomy_utax.txt > otus_tax_blast.fna
~~~
OTU table construction
~~~
usearch81_64_new -usearch_global seqs_dot.fna -db otus_tax.fna -strand plus -id 0.985 -otutabout otu_table.txt
~~~
>	868329 / 933570 mapped to OTUs (93.0%) 

## BIOM table construction

Convert UTAX taxonomy to phyloseq compatible
- Separate taxonomy and sequence counts
~~~
cut -f 63- otu_table.txt > tax.txt
cut -f 1-62 otu_table.txt > abund.txt
~~~

- Delete taxonomy level letters, change field separator from "," to ";", delete fileds containing "uncultured*", fill empty fields with "u_(lowest assigned level)".
~~~
sed -E 's/d://g; s/p://g; s/c://g; s/o://g; s/f://g; s/g://g; s/s://g' tax.txt| awk -F"," '{ if ($2 == "") print $1",u_"$1;  else print $0}' | awk -F"," '{ if ($3 == "") print $1","$2",u_"$2;  else print $0}' | awk -F"," '{ if ($4 == "") print $1","$2","$3",u_"$3;  else print $0}' | awk -F"," '{ if ($5 == "") print $1","$2","$3","$4",u_"$4;  else print $0}' | awk -F"," '{ if ($6 == "") print $1","$2","$3","$4","$5",u_"$5;  else print $0}' | awk -F"," '{ if ($7 == "") print $1","$2","$3","$4","$5","$6",u_"$6;  else print $0}' | sed -E '1s/^.*$/taxonomy/; s/,/;/g' | sed 's/\(u_\)\1\{1,\}/u_/g' > tax_phyloseq.txt
~~~

- Merge sequence counts and phyloseq comaptible taxonomy
~~~
paste abund.txt tax_phyloseq.txt > otu_table_phyloseq.txt
~~~

Discard non-fungal OTUs
~~~
sed -i '3,${/Fungi/!d;}' otu_table.txt
~~~

Create and summarize biom file
~~~
biom convert -i otu_table_phyloseq.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o otu_table_phyloseq.biom

biom summarize-table -i otu_table_phyloseq.biom -o otu_table_phyloseq_summary.txt
~~~

> Num samples: 61
Num observations: 3098
Total count: 922919
Table density (fraction of non-zero values): 0.113

> Counts/sample summary:
 Min: 13.0
 Max: 38220.0
 Median: 13134.000
 Mean: 15129.820
 Std. dev.: 9512.690
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

> Counts/sample detail:
577: 13.0
529: 1431.0
570: 1621.0
533: 1786.0
554: 2995.0
581: 3668.0
561: 3979.0
541: 5150.0
553: 5611.0
537: 5744.0
530: 5979.0
555: 6190.0
578: 6679.0
521: 7091.0
573: 7196.0
580: 7362.0
568: 7476.0
545: 7630.0
522: 8558.0
557: 8689.0
572: 8784.0
525: 8785.0
546: 8944.0
579: 9368.0
582: 9473.0
524: 9560.0
562: 9881.0
528: 10117.0
556: 12343.0
566: 12740.0
523: 13134.0
575: 15859.0
536: 16197.0
527: 16572.0
574: 16895.0
532: 17009.0
526: 17341.0
534: 17894.0
551: 17952.0
531: 19066.0
571: 19998.0
565: 20180.0
543: 20404.0
550: 21395.0
560: 21716.0
567: 22612.0
539: 22775.0
523: 13134.0
575: 15859.0
536: 16197.0
527: 16572.0
574: 16895.0
532: 17009.0
526: 17341.0
534: 17894.0
551: 17952.0
531: 19066.0
571: 19998.0
565: 20180.0
543: 20404.0
550: 21395.0
560: 21716.0
567: 22612.0
539: 22775.0
540: 23566.0
538: 23587.0
548: 25303.0
547: 25345.0
576: 26527.0
535: 26788.0
559: 26994.0
544: 27067.0
563: 28675.0
552: 29141.0
558: 29267.0
549: 33308.0
564: 37289.0
542: 38220.0


### SOFTWARE USED
- QIIME 1.9.1
- PRINSEQ-lite 0.20.4
- sed (GNU sed) 4.2.2
- GNU Awk 4.1.1, API: 1.1 (GNU MPFR 3.1.2-p3, GNU MP 6.0.0)
- cut (GNU coreutils) 8.23
- paste (GNU coreutils) 8.23
- biom, version 2.1.5
- ITSx version 1.0.11
