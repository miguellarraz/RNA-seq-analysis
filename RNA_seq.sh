
#!/bin/bash
#Directories I will be needing
echo "Starting analysis"
mkdir fastqc_results && mkdir index_files && mkdir aligned_reads && mkdir counts && mkdir bamfiles && mkdir final_counts && mkdir fastqc_reports && mkdir groupwise_comparisons
#to align the reads, I need to make an index of the reference genome in a directory called index_files, only need to do this once, so outside$
bowtie2-build -q --threads 60 /localdisk/home/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz index_files/index_ref_file && echo "Created index files"
for file in /localdisk/home/data/BPSM/AY21/fastq/*.fq.gz #looping for every file, and running fastqc on 60 threads
do
#-f flag used, to specify the input format, in this case fastq
fastqc -t 60 -f fastq -q $file -o $PWD/fastqc_results 
done
echo "Fastqc analysis completed, creating report..."
#making a file that has both forwards and reverse reads on the same line, so that I can use a single variable for both. I do this by first writing down all the forwards reads(they end in 1), and then appending the reverse reads to the same rows as the corresponding forwards read
touch forward_and_reverse_reads.txt && touch forward_reads.txt && touch reverse_reads.txt
for file in /localdisk/home/data/BPSM/AY21/fastq/*1.fq.gz
do
echo "${file:37}" >> forward_reads.txt
done
for file in /localdisk/home/data/BPSM/AY21/fastq/*2.fq.gz
do
echo "${file:37}" >> reverse_reads.txt
done
paste -d "" forward_reads.txt reverse_reads.txt > forward_and_reverse_reads.txt
#Making the FASTQC report
counter_0=0 #two counters to help with the file name
counter_1=1
#First, add the headers for the fastqc report
echo "Sample" >> fastqc_categories.txt && echo "Basic Statistics" >> fastqc_categories.txt && echo "Per base sequence quality" >> fastqc_categories.txt && echo "Per sequence quality scores" >> fastqc_categories.txt && echo "Per base sequence content" >> fastqc_categories.txt && echo "Per sequence GC content" >> fastqc_categories.txt && echo "Per base N content" >> fastqc_categories.txt && echo "Sequence Length Distribution" >> fastqc_categories.txt && echo "Sequence Duplication Levels" >> fastqc_categories.txt && echo "Overrepresented sequences" >> fastqc_categories.txt && echo "Adapter Content" >> fastqc_categories.txt
cp fastqc_categories.txt fastqc_reports/fastqc_report_0.txt #copying this file into the first version of the report, so that I can append the scores one by one
cat forward_reads.txt reverse_reads.txt > all_reads.txt #A file that contains all the reads, easier to loop
cd fastqc_results
unzip -qq '*.zip' #unzip all the files, easier if I am in the directory with all the files
cd ..
cat all_reads.txt | while read filename
do
current_file=$(echo ${filename:0:15}) #The name of the current file
echo "$current_file" > temp_fastqc_file.txt #as a header
cut -f 1 fastqc_results/"$current_file"_fastqc/summary.txt >> temp_fastqc_file.txt # writing the summary into this temporary file
paste -d'\t' fastqc_reports/fastqc_report_$counter_0.txt temp_fastqc_file.txt > fastqc_reports/fastqc_report_$counter_1.txt # Doing this so that I can sequentially append a summary to make the full report
rm -fr temp_fastqc_file.txt #Don't need this anymore
((counter_0=counter_0+1))
((counter_1=counter_1+1))
done
finished_report=$(wc -l  all_reads.txt | cut -d ' ' -f1) #This is so that I knoe how many files I made, so that I can chose the last one, which should be the full report
mv fastqc_reports/*$finished_report.txt fastqc_report.txt #full report
rm -fr fastqc_reports && rm -fr all_reads && rm -fr fastqc_categories && rm -fr forward_reads && rm -fr reverse_reads && rm -fr forward_and_reverse_reads
echo "Report finished. Aligning reads and producing gene counts..."
#Making a file that shows all the gene counts
#make the first gene counts file, this only has the gene names and descriptions
touch genes.txt
echo "Gene	Gene_description">>genes.txt
cut -f 4,5 /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed >>genes.txt
#make a copy for the final gene counts, this is the format I will be using later(gene_counts_0.txt)
cp genes.txt final_counts/gene_counts_0.txt
cat forward_and_reverse_reads.txt | while read line #Using the file I created before with the forwards and reverse reads in the same line for each culture
do
input_file_1=$(echo ${line:0:21}) #in this file, the first part of each line is the forwards read, and the second the reverse read
input_file_2=$(echo ${line:21})
sam_name=${line:0:13} #This is the sample name, which I will use to name the sam file
#the alignment itself, with output in aligned_reads
bowtie2 -p 60 -x index_files/index_ref_file -1 /localdisk/home/data/BPSM/AY21/fastq/$input_file_1 -2 /localdisk/home/data/BPSM/AY21/fastq/$input_file_2 -S aligned_reads/$sam_name.aligned.sam 2> stats.file
echo "$sam_name" >> sam_filenames.txt
lastline=$(tail -1 stats.file) #Here I output the info on the allignment, so that I can pich the % allignment and output a message to a log file if the value is smaller than a threshold
success=$(echo ${lastline:0:1}) #first character of the percentage
alignmentrate=$(echo ${lastline:0:6})
if [ "$success" -lt "8" ] #threshold. I choose 80%
then
   echo "Culture "$sam_name" has a rate of alignment of only "$alignmentrate". You might want to test for contamination using BLAST" >> logfile.log
fi 
done
echo "Alignment finished"
i=0
p=1
#convert output to indexed bam
sort -t$'\t' -k1 genes.txt > genes.sorted #I want to sort the genes.txt file to make sure all reads are assigned to the correct gene
cat sam_filenames.txt | while read samfile #sam_filenames comes from part of each line in forwards_and_reverse_reads, and it is the name of each replicate.
do
samtools view -b aligned_reads/$samfile.aligned.sam  > bamfiles/$samfile.aligned.bam && samtools sort bamfiles/$samfile.aligned.bam > bamfiles/$samfile.aligned.sorted.bam && samtools index bamfiles/$samfile.aligned.sorted.bam #Making bam dile with samtools, then sorting it, and indexing it
bedtools bamtobed -i bamfiles/$samfile.aligned.sorted.bam > bamfiles/$samfile.aligned.sorted.bed #making bed
#Using bedtools coverage with the counts flag to produce counts for every gene
bedtools coverage -counts -a /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b bamfiles/$samfile.aligned.sorted.bed > counts/$samfile.counts.txt
#Now I want to merge the gene count file, but need to sort, as you need it before the merging, thats why I sorted genes.txt before
sort -t$'\t' -k4 counts/$samfile.counts.txt > counts/$samfile.counts.sorted.txt
#cut the gene id and the count columns
cut -f 4,6 counts/$samfile.counts.sorted.txt > temp_count_file.txt
sed -i "1i Gene""	$samfile" temp_count_file.txt #Adding headers
join -t$'\t' -e 'NaN' genes.sorted temp_count_file.txt > temp_gene_counts.txt #Joining both files, probably not necessary, but I want to make sure there are the same number of genes after bedtools coverage. -e flag is in case a gene is not present 
#quick error message to the logfile if the output of bedtools coverage is missing genes
genelines1=$(wc -l  genes.txt | cut -d ' ' -f1)
genelines2=$(wc -l  temp_gene_counts.txt | cut -d ' ' -f1)
if [ $genelines1 != $genelines2 ]
then
   echo "Error in replicate "$samfile"" >> logfile.log

fi 
cut -f 3 temp_gene_counts.txt > temp_gene_counts_single_columns.txt #single column containing the counts of a replicate
paste -d'\t' final_counts/gene_counts_$i.txt temp_gene_counts_single_columns.txt > final_counts/gene_counts_$p.txt #pasting the counts onto the previous counts file. First one contains only the genes and descriptions
((i=i+1))
((p=p+1))
done
#Select the last file, which is the one containing all the gene counts
last_file=$(wc -l  sam_filenames.txt | cut -d ' ' -f1)
mv final_counts/*$last_file.txt replicate_gene_expressions.txt
echo "Producing mean counts..."
#grouping into replicate groups
mkdir mean_counts
tail -n+2 /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | awk '{print $2}'| sort | uniq > mean_counts/samples.txt #samples file, so that more samples can be added
cat mean_counts/samples.txt | while read sample #for each sample
do
# 0 hours uninduced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | awk -F '\t' '$4 == "0" {print $1}' > mean_counts/$sample.0hr.replicates.txt #all the names of the replicates of a single sample, zero hours
   cat mean_counts/$sample.0hr.replicates.txt | while read replicate #Read the file that contains all the zero hour cultures for that sample
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1) #get the column position for this replicate
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp #make a file with this column
   done
   paste mean_counts/*temp > mean_counts/$sample.0hr.uninduced.samples #paste all columns into the same file
   rm -fr mean_counts/*temp && rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.0hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.0hr.uninduced.mean #the mean of each column, and output it as an extra column in the file
   sed -i "1i""$sample"".0hr.uninduced" mean_counts/$sample.0hr.uninduced.mean #and add column title
   rm -fr mean_counts/*samples
#24 hours induced, same as before, comments awould be the same as above
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Induced"| awk -F '\t' '$4 == "24" {print $1}' > mean_counts/$sample.24hr.induced.replicates.txt
   cat mean_counts/$sample.24hr.induced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.24hr.induced.samples
   rm -fr mean_counts/*temp && rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.24hr.induced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.24hr.induced.mean
   sed -i "1i""$sample"".24hr.induced" mean_counts/$sample.24hr.induced.mean
   rm -fr mean_counts/*samples
#24 hours uninduced, same as before
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Uninduced"| awk -F '\t' '$4 == "24" {print $1}' > mean_counts/$sample.24hr.uninduced.replicates.txt
   cat mean_counts/$sample.24hr.uninduced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.24hr.uninduced.samples
   rm -fr mean_counts/*temp && rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.24hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.24hr.uninduced.mean
   sed -i "1i""$sample"".24hr.uninduced" mean_counts/$sample.24hr.uninduced.mean
   rm -fr mean_counts/*samples
#48 hours induced, same as before
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Induced"| awk -F '\t' '$4 == "48" {print $1}' > mean_counts/$sample.48hr.induced.replicates.txt
   cat mean_counts/$sample.48hr.induced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.48hr.induced.samples
   rm -fr mean_counts/*temp && rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.48hr.induced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.48hr.induced.mean
   sed -i "1i""$sample"".48hr.induced" mean_counts/$sample.48hr.induced.mean
   rm -fr mean_counts/*samples
#48 hours uninduced, same as before
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Uninduced"| awk -F '\t' '$4 == "48" {print $1}' > mean_counts/$sample.48hr.uninduced.replicates.txt
   cat mean_counts/$sample.48hr.uninduced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.48hr.uninduced.samples
   rm -fr mean_counts/*temp && rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.48hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.48hr.uninduced.mean
   sed -i "1i""$sample"".48hr.uninduced" mean_counts/$sample.48hr.uninduced.mean
   rm -fr mean_counts/*samples
done
paste mean_counts/*mean > mean_counts/replicate_mean_counts.txt #joining all means together, and adding gene names below
paste genes.txt mean_counts/replicate_mean_counts.txt > replicate_mean_counts.txt
cp mean_counts/samples.txt samples.txt
echo "Mean counts created, calculating fold change for all pairwise comparisons"
#Fold change
#calculating all possible comparisons, and writing them down to a file
head -1 replicate_mean_counts.txt | tr "\t" "\n" | grep -v Gene | grep -v Gene_description > all_groups.txt
cat all_groups.txt | while read i
do
  cat all_groups.txt | while read j
  do
    if [ "$i" \< "$j" ] #find all combinations for each file
    then
     echo "$i:$j" >> all_combinations.txt
     echo "$j:$i" >> all_combinations.txt
    fi
  done
done
rm -fr all_groups.txt
#Now to find the fold change for every possible comparison
cat all_combinations.txt | while read STR
do
   comp1=$(echo $STR | cut -f1 -d ':') #use the separator : to distinguish between files
   comp2=$(echo $STR | cut -f2 -d ':')

   column1=$(head -1 replicate_mean_counts.txt | tr '\t' '\n' | cat -n | grep $comp1 | cut -d$'\t' -f1) #find the positions of the columns, cutting and pasting them to a new file
   column2=$(head -1 replicate_mean_counts.txt | tr '\t' '\n' | cat -n | grep $comp2 | cut -d$'\t' -f1)
   cut -f $column1 replicate_mean_counts.txt > groupwise_comparisons/mean.counts.$comp1 && cut -f $column2 replicate_mean_counts.txt > groupwise_comparisons/mean.counts.$comp2 #cut the columns for both comparisons
   paste -d'\t' groupwise_comparisons/mean.counts.$comp1 groupwise_comparisons/mean.counts.$comp2 > groupwise_comparisons/comparison_columns #paste them together
   sed -i '1s/$/  fold_change/' groupwise_comparisons/comparison_columns #add a header to the next column
   awk -v OFS='\t' 'NR!=1 {$3 = ($2 != 0) ? sprintf("%.3f", $1/$2) : "NAN"}1' groupwise_comparisons/comparison_columns > groupwise_comparisons/$comp1.vs.$comp2.single_column #Calculating fold change as y/x, and adding NAN when the initial value is zero
   paste -d'\t' genes.sorted groupwise_comparisons/$comp1.vs.$comp2.single_column | tail -n+2 | sort -t$'\t' -n -k5 -r > groupwise_comparisons/$comp1.vs.$comp2 #adding gene descriptions and names
   sed -i "1i""Gene  Gene_description  ""$comp1""  ""$comp2""  ""fold_change" groupwise_comparisons/$comp1.vs.$comp2 #Adding headers
   rm -fr groupwise_comparisons/mean* && rm -fr groupwise_comparisons/comparison* && rm -fr groupwise_comparisons/*single_column
done
echo "Done!!"
#Remove messy files and directories
rm -fr forward_and_reverse_reads.txt && rm -fr forward_reads.txt && rm -fr reverse_reads.txt && rm -fr temp_gene_counts_single_columns.txt && rm -fr temp_gene_counts.txt && rm -fr temp_count_file.txt && rm -fr genes.txt && rm -fr genes.sorted.txt && rm -fr sam_filenames.txt && rm -fr all_combinations.txt
rm -fr counts && rm -fr final_counts && rm -fr mean_counts && rm -fr stats.file && rm -fr all_reads.txt && rm -fr fastqc_categories.txt && rm -fr samples.txt && rm -fr genes.sorted
