#Individual library cleaning steps
#Acer VanWallendael adapted from https://github.com/Gian77/Amplicon-Pipelines-hpcc

#unzip reads
gunzip -c Lib14_phyllos_S1_L001_R1_001.fastq.gz > r1.fq
gunzip -c Lib14_phyllos_S1_L001_I1_001.fastq.gz  >  i1.fq
gunzip -c Lib14_phyllos_S1_L001_R2_001.fastq.gz > r2.fq

mkdir phyllos_all
cp *.fq ./phyllos_all
cd phyllos_all

#count the reads
for i in *.fq
do echo "$i : `grep -c "^+$" $i`"
done > reads_raw.txt

#load qiime
conda create -n qiime1 qiime mock nose -c bioconda

source activate qiime1
print_qiime_config.py -t

#get some QC
fastqc ./r1.fq
fastqc ./r2.fq

#need to set up the map file
#from drive folder
#mv ../Barcodes_phyllos.tsv map_ITS.txt

#demultiplex
validate_mapping_file.py -m map_ITS_corrected.txt -o mapping_output
#needs # before sampleid
#no _ in sample IDs
#can use corrected
cp ./mapping_output/map_ITS_corrected.txt map_ITS1.txt

validate_mapping_file.py -m map_ITS1.txt -o mapping_output1

split_libraries_fastq.py -i r1.fq -m map_ITS_corrected.txt -b i1.fq --store_demultiplexed_fastq --barcode_type 10 -o demultiplexed_R1/ --rev_comp_mapping_barcodes -q 19 --max_barcode_errors 0
split_libraries_fastq.py -i r2.fq -m map_ITS_corrected.txt -b i1.fq --store_demultiplexed_fastq --barcode_type 10 -o demultiplexed_R2/ --rev_comp_mapping_barcodes -q 19 --max_barcode_errors 0

#clean for fast*
sed '/^>/s/_/\./g' "demultiplexed_R1/seqs.fna" > "demultiplexed_R1/seqs_new.fasta"
sed '/^>/s/_/\./g' "demultiplexed_R2/seqs.fna" > "demultiplexed_R2/seqs_new.fasta"
sed '/^@/s/_/\./g' "demultiplexed_R1/seqs.fastq" > "demultiplexed_R1/seqs_new.fastq"
sed '/^@/s/_/\./g' "demultiplexed_R2/seqs.fastq" > "demultiplexed_R2/seqs_new.fastq"

#split into samples
split_sequence_file_on_sample_ids.py -i demultiplexed_R1/seqs.fastq --file_type fastq -o splitted/
for f in splitted/*fastq ; do NEW=${f%.fastq}_R1.fastq; mv ${f} "${NEW}"; done
split_sequence_file_on_sample_ids.py -i demultiplexed_R2/seqs.fastq --file_type fastq -o splitted_R2/
for f in splitted_R2/*fastq ; do NEW=${f%.fastq}_R2.fastq; mv ${f} "${NEW}"; done

mv splitted_R2/*R2.fastq splitted/

mkdir stripped/
#conda create -n cutadapt
#then restart connection
#conda activate cutadaptenv
#python3 -m pip install --user --upgrade cutadapt

module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load cutadapt/2.1-Python-3.6.6

cutadapt -g CTTGGTCATTTAGAGGAAGTAA -e 0.01 --discard-untrimmed --match-read-wildcards demultiplexed_R1/seqs_new.fastq > stripped/stripped_R1.fastq
cutadapt -g TCCTCCGCTTATTGATATGC -e 0.01 --discard-untrimmed --match-read-wildcards demultiplexed_R2/seqs_new.fastq > stripped/stripped_R2.fastq

#conda deactivate
module load GCC/5.4.0-2.26  OpenMPI/1.10.3-CUDA
module load seqtk

seqtk sample -s100 stripped/stripped_R1.fastq 500 > stripped/sub_stripped_R1.fastq
seqtk seq -aQ64 stripped/sub_stripped_R1.fastq > stripped/sub_stripped_R1.fasta

seqtk sample -s100 stripped/stripped_R2.fastq 500 > stripped/sub_stripped_R2.fastq  	 
seqtk seq -aQ64 stripped/sub_stripped_R2.fastq > stripped/sub_stripped_R2.fasta

mkdir stats

for fastq in stripped/*.fastq; do echo "$fastq : `grep -c "^+$" $fastq`"; done > stats/stripped.counts

module load FastQC/0.11.7-Java-1.8.0_162
module load vsearch/2.9.1

mkdir filtered

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter stripped/stripped_R1.fastq -fastq_stripleft 45 -fastqout stripped/trimmed_R1.fastq
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter stripped/stripped_R2.fastq -fastq_stripleft 36 -fastqout stripped/trimmed_R2.fastq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter stripped/trimmed_R1.fastq -fastq_maxee 1.0 -fastq_trunclen 200 -fastq_maxns 0 -fastqout filtered/filtered_R1.fastq
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter stripped/trimmed_R2.fastq -fastq_maxee 1.0 -fastq_trunclen 200 -fastq_maxns 0 -fastqout filtered/filtered_R2.fastq
fastqc ./filtered/filtered_R1.fastq
fastqc ./filtered/filtered_R2.fastq