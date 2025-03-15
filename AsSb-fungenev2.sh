#!/usr/bin/env sh
#正常直接运行的版本
INPUT=$1

mkdir ${INPUT}_outputs
mkdir ${INPUT}_tmp
# Data cleaning using fastp
fastp -i ${INPUT}_1.fastq -I ${INPUT}_2.fastq -q 30 -l 25 -o ${INPUT}_P1.fastq -O ${INPUT}_P2.fastq -w 30 -h ${INPUT}_fastp.html -j ${INPUT}_fastp.json --detect_adapter_for_pe --trim_poly_g --trim_poly_x 

rm ${INPUT}_1.fastq
rm ${INPUT}_2.fastq
# Assembly using megahit
megahit -1 ${INPUT}_P1.fastq -2 ${INPUT}_P2.fastq -o ${INPUT}_assembly -t 12 -m 0.9 --k-step 10
# --k-min 31 --k-max 141 

rm ${INPUT}_P1.fastq
rm ${INPUT}_P2.fastq

bowtie2-build ${INPUT}_assembly/final.contigs.fa ${INPUT}_tmp/${INPUT}.index
bowtie2 -q ${INPUT}_cleaned.fastq -x ${INPUT}_tmp/${INPUT}.index -p 6 | samtools view -Sb | samtools sort > ${INPUT}_tmp/${INPUT}_automapped_sorted_bam

prodigal -i ${INPUT}_assembly/final.contigs.fa -o ${INPUT}_tmp/${INPUT}_genes.gff -f gff -a ${INPUT}_tmp/${INPUT}_proteins.faa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${INPUT}_tmp/${INPUT}_proteins.faa > ${INPUT}_tmp/${INPUT}_proteins2.faa
awk '{print $NF}' ${INPUT}_tmp/${INPUT}_proteins2.faa > ${INPUT}_tmp/${INPUT}_protseq.txt
awk -F '#' '{print $1 }' ${INPUT}_tmp/${INPUT}_proteins2.faa  > ${INPUT}_tmp/${INPUT}_geneid.txt
paste ${INPUT}_tmp/${INPUT}_geneid.txt ${INPUT}_tmp/${INPUT}_protseq.txt | cut -f 1,2 > ${INPUT}_tmp/${INPUT}_prot1.txt
sed 's/>//' ${INPUT}_tmp/${INPUT}_prot1.txt > ${INPUT}_tmp/${INPUT}_prot2.txt
awk '{sub("-", "", $2); print}' < ${INPUT}_tmp/${INPUT}_prot2.txt > ${INPUT}_tmp/${INPUT}_prot2b.txt
awk '{gsub(/*$/,""); print}' ${INPUT}_tmp/${INPUT}_prot2b.txt > ${INPUT}_tmp/${INPUT}_prot2c.txt
echo -e "gene_id prot" | cat - ${INPUT}_tmp/${INPUT}_prot2c.txt  > ${INPUT}_tmp/${INPUT}_prot3.txt

featureCounts -t CDS -o ${INPUT}_tmp/${INPUT}_counts.tsv -g ID -a ${INPUT}_tmp/${INPUT}_genes.gff ${INPUT}_tmp/${INPUT}_automapped_sorted_bam
sed -e '1,2d' ${INPUT}_tmp/${INPUT}_counts.tsv > ${INPUT}_tmp/${INPUT}_counts.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $2}' > ${INPUT}_tmp/${INPUT}_id1.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $1}' | cut -d'_' -f2 > ${INPUT}_tmp/${INPUT}_id2.txt
paste -d":" ${INPUT}_tmp/${INPUT}_id1.txt ${INPUT}_tmp/${INPUT}_id2.txt | sed 's/:/_/g' > ${INPUT}_tmp/${INPUT}_id3.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $7}' > ${INPUT}_tmp/${INPUT}_counts1.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $6}' > ${INPUT}_tmp/${INPUT}_length.txt
paste ${INPUT}_tmp/${INPUT}_id3.txt ${INPUT}_tmp/${INPUT}_counts1.txt ${INPUT}_tmp/${INPUT}_length.txt > ${INPUT}_tmp/${INPUT}_counts2.txt
echo -e "gene_id seq length" | cat - ${INPUT}_tmp/${INPUT}_counts2.txt  > ${INPUT}_tmp/${INPUT}_counts3.txt

# searching arsA ATPase homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsA.txt --tblout ${INPUT}_outputs/${INPUT}_arsA_hmmer.out hmmdb/arsGenes/ArsA_ATPase-PF02374.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsA_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsA_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsA_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsA_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsA_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsA_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsA_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsA_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsA_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsA_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsA_counts.txt ${INPUT}_tmp/${INPUT}_arsA_seq.txt | sed 's/$/ arsA_hom/'  > ${INPUT}_outputs/${INPUT}_arsA_homologs.txt

# searching arsB homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsB.txt --tblout ${INPUT}_outputs/${INPUT}_arsB_hmmer.out hmmdb/arsGenes/ArsB-PF02040.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsB_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsB_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsB_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsB_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsB_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsB_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsB_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsB_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsB_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsB_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsB_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsB_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsB_counts.txt ${INPUT}_tmp/${INPUT}_arsB_seq.txt | sed 's/$/ arsB_hom/'  > ${INPUT}_outputs/${INPUT}_arsB_homologs.txt

# searching arsC homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsC.txt --tblout ${INPUT}_outputs/${INPUT}_arsC_hmmer.out hmmdb/arsGenes/ArsC-PF03960.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsC_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsC_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsC_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsC_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsC_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsC_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsC_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsC_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsC_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsC_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsC_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsC_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsC_counts.txt ${INPUT}_tmp/${INPUT}_arsC_seq.txt | sed 's/$/ arsC_hom/'  > ${INPUT}_outputs/${INPUT}_arsC_homologs.txt

# searching arsD homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsD.txt --tblout ${INPUT}_outputs/${INPUT}_arsD_hmmer.out hmmdb/arsGenes/ArsD-PF06953.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsD_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsD_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsD_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsD_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsD_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsD_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsD_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsD_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsD_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsD_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsD_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsD_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsD_counts.txt ${INPUT}_tmp/${INPUT}_arsD_seq.txt | sed 's/$/ arsD_hom/'  > ${INPUT}_outputs/${INPUT}_arsD_homologs.txt

# searching arsR homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsR.txt --tblout ${INPUT}_outputs/${INPUT}_arsR_hmmer.out hmmdb/arsGenes/ArsR-PF09824.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsR_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsR_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsR_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsR_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsR_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsR_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsR_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsR_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsR_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsR_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsR_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsR_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsR_counts.txt ${INPUT}_tmp/${INPUT}_arsR_seq.txt | sed 's/$/ arsR_hom/'  > ${INPUT}_outputs/${INPUT}_arsR_homologs.txt

# searching arsA hsp20 homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_ArsAhsp.txt --tblout ${INPUT}_outputs/${INPUT}_ArsAhsp_hmmer.out hmmdb/arsGenes/ArsA_HSP20-PF17886.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_ArsAhsp_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_ArsAhsp_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_ArsAhsp_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins2.faa > ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_ArsAhsp_proteins3.faa > ${INPUT}_tmp/${INPUT}_ArsAhsp_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_ArsAhsp_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_ArsAhsp_counts.txt
paste ${INPUT}_tmp/${INPUT}_ArsAhsp_counts.txt ${INPUT}_tmp/${INPUT}_ArsAhsp_seq.txt | sed 's/$/ ArsAhsp_hom/'  > ${INPUT}_outputs/${INPUT}_ArsAhsp_homologs.txt

# searching arsP homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsP.txt --tblout ${INPUT}_outputs/${INPUT}_arsP_hmmer.out hmmdb/arsGenes/ArsP-PF11449.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsP_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsP_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsP_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsP_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsP_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsP_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsP_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsP_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsP_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsP_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsP_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsP_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsP_counts.txt ${INPUT}_tmp/${INPUT}_arsP_seq.txt | sed 's/$/ arsP_hom/'  > ${INPUT}_outputs/${INPUT}_arsP_homologs.txt

# searching arsP_1 homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsP_1.txt --tblout ${INPUT}_outputs/${INPUT}_arsP_1_hmmer.out hmmdb/arsGenes/ArsP_1-PF03773.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsP_1_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsP_1_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsP_1_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsP_1_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsP_1_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsP_1_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsP_1_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsP_1_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsP_1_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsP_1_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsP_1_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsP_1_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsP_1_counts.txt ${INPUT}_tmp/${INPUT}_arsP_1_seq.txt | sed 's/$/ arsP_1_hom/'  > ${INPUT}_outputs/${INPUT}_arsP_1_homologs.txt


# searching ars2 homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_ars2.txt --tblout ${INPUT}_outputs/${INPUT}_ars2_hmmer.out hmmdb/arsGenes/Ars2-PF04959.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_ars2_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_ars2_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_ars2_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_ars2_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_ars2_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_ars2_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_ars2_proteins2.faa > ${INPUT}_tmp/${INPUT}_ars2_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_ars2_proteins3.faa > ${INPUT}_tmp/${INPUT}_ars2_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_ars2_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_ars2_counts.txt
paste ${INPUT}_tmp/${INPUT}_ars2_counts.txt ${INPUT}_tmp/${INPUT}_ars2_seq.txt | sed 's/$/ ars2_hom/'  > ${INPUT}_outputs/${INPUT}_ars2_homologs.txt


# searching ars2N homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_ars2N.txt --tblout ${INPUT}_outputs/${INPUT}_ars2N_hmmer.out hmmdb/arsGenes/Ars2N-PF12066.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_ars2N_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_ars2N_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_ars2N_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_ars2N_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_ars2N_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_ars2N_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_ars2N_proteins2.faa > ${INPUT}_tmp/${INPUT}_ars2N_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_ars2N_proteins3.faa > ${INPUT}_tmp/${INPUT}_ars2N_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_ars2N_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_ars2N_counts.txt
paste ${INPUT}_tmp/${INPUT}_ars2N_counts.txt ${INPUT}_tmp/${INPUT}_ars2N_seq.txt | sed 's/$/ ars2N_hom/'  > ${INPUT}_outputs/${INPUT}_ars2N_homologs.txt

#Consider add--cut_tc for following code
# searching arrA homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arrA.txt --tblout ${INPUT}_outputs/${INPUT}_arrA_hmmer.out hmmdb/arrA.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arrA_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arrA_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arrA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arrA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arrA_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arrA_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arrA_proteins2.faa > ${INPUT}_tmp/${INPUT}_arrA_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arrA_proteins3.faa > ${INPUT}_tmp/${INPUT}_arrA_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arrA_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arrA_counts.txt
paste ${INPUT}_tmp/${INPUT}_arrA_counts.txt ${INPUT}_tmp/${INPUT}_arrA_seq.txt | sed 's/$/ arrA_hom/'  > ${INPUT}_outputs/${INPUT}_arrA_homologs.txt

# searching arsH homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arsH.txt --tblout ${INPUT}_outputs/${INPUT}_arsH_hmmer.out hmmdb/arsGenes/arsH.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arsH_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arsH_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arsH_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arsH_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arsH_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arsH_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arsH_proteins2.faa > ${INPUT}_tmp/${INPUT}_arsH_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arsH_proteins3.faa > ${INPUT}_tmp/${INPUT}_arsH_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arsH_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arsH_counts.txt
paste ${INPUT}_tmp/${INPUT}_arsH_counts.txt ${INPUT}_tmp/${INPUT}_arsH_seq.txt | sed 's/$/ arsH_hom/'  > ${INPUT}_outputs/${INPUT}_arsH_homologs.txt

# searching arrA homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_arrA.txt --tblout ${INPUT}_outputs/${INPUT}_arrA_hmmer.out hmmdb/arrA.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_arrA_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_arrA_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_arrA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_arrA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_arrA_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_arrA_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_arrA_proteins2.faa > ${INPUT}_tmp/${INPUT}_arrA_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_arrA_proteins3.faa > ${INPUT}_tmp/${INPUT}_arrA_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_arrA_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_arrA_counts.txt
paste ${INPUT}_tmp/${INPUT}_arrA_counts.txt ${INPUT}_tmp/${INPUT}_arrA_seq.txt | sed 's/$/ arrA_hom/'  > ${INPUT}_outputs/${INPUT}_arrA_homologs.txt




# searching glpF homologs
hmmsearch -o ${INPUT}_tmp/${INPUT}_glpF.txt --tblout ${INPUT}_outputs/${INPUT}_glpF_hmmer.out hmmdb/glpF.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_glpF_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_glpF_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_glpF_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_glpF_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_glpF_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_glpF_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_glpF_proteins2.faa > ${INPUT}_tmp/${INPUT}_glpF_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_glpF_proteins3.faa > ${INPUT}_tmp/${INPUT}_glpF_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_glpF_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_glpF_counts.txt
paste ${INPUT}_tmp/${INPUT}_glpF_counts.txt ${INPUT}_tmp/${INPUT}_glpF_seq.txt | sed 's/$/ glpF_hom/'  > ${INPUT}_outputs/${INPUT}_glpF_homologs.txt


# Cleaning folder from temporary files
# rm -r ${INPUT}_tmp
