#A bash script for generating vial analytic database
#Input file must be a fasta nuclotide alignment aligned correctly to the GFF coordinates

#clean fasta headers to remove any un-allowed characters

cat aligned_input.fasta | tr -d '/' | tr -d '|' | tr -d ' ' | tr '-' 'N' > ./split_fasta/clean_input

cd split_fasta

#split tha multifata file in to single files
faidx -x ./clean_input

rm clean_input clean_input.fai

#generate gff files for each single fasta file
for i in `ls` ; do cat ../SARS_CoV_2.gff | sed "s/MT188341/$i/g"  > ../gff/$i.gff  ; done

#extract regions listed in the GFF file

for i in `ls` ; do bedtools getfasta -fi ./$i -bed ../gff/$i.gff -fo ./../sequence_processing/$i -name ; done

cd ../sequence_processing

#extract the individaul protein names
for i in `ls` ; do cat $i ; done > all_protein_combined.fasta

cat all_protein_combined.fasta | grep \> | sort -u | sed 's/>//g' > ../protein_names.txt

#collate individual protein data in to a single file
while read -r line; do fgrep -w  $line  all_protein_combined.fasta -A 1  > ../combined_protein_cds/$line.fasta ; done < ../protein_names.txt

cd ../combined_protein_cds/

#translate sequences
for i in `ls *.fasta`; do transeq -trim -sequence $i -outseq ../translated/$i ; done

#generate cannonical sequences
cd ../translated

#for i in `ls *.fasta` ; do cons -sequence $i -outseq ../cannonical_sequences/$i -name $i ; done

#rename the protein fasta files
cat ../aligned_input.fasta | grep \> > ./sequence_names.txt

#trim fasta names off

#replacement loop
for i in `ls *.fasta` ; do awk 'NR==FNR{names[NR]=$0; next} /^>/{$1=""names[++c]}1' sequence_names.txt $i > ./$i.renamed; done

rm sequence_names.txt *.fasta

rename s/\orf1ab_// *.renamed

rename s/.fasta.renamed// *.renamed

rename s/\>// *

for i in `ls` ; do pyfasta split -n 2 $i ; done

for i in `ls -I "*.*"` ; do mafft --auto --thread -1  --keeplength --addfragments $i.0 ../wuhan_reference_protein_seq/$i > $i.aligned ; done

for i in `ls -I "*.*"` ; do mafft --auto --thread -1  --keeplength --addfragments $i.1 ../wuhan_reference_protein_seq/$i > $i.1.aligned ; done

rm *.flat *.gdx

for i in `ls -I "*.*"` ; do cat $i.aligned $i.1.aligned > $i ; done

rm *.*

for i in `ls` ; do python ../NON_SYN_CALL.py $i 8 ; done

cp *.csv ../mutations

cd ../mutations

rename s/.AminoAcidCounts.csv// *.csv

sed -i 's/North America/North-America/g' *

sed -i 's/South America/South-America/g' *

#for i in `ls` ; do cat $i | tr ',' '\t' | awk -v name=$i  '{print name"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8(.*)}' > $i.genecolumn; done
for i in `ls` ; do cat $i | tr ',' '\t' | sed 's/\,\EPI/EPI/g'  |  awk -F ","  -v name=$i  '{print name"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > $i.genecolumn ; done

rename -f s/.genecolumn// *.genecolumn

for i in `ls`; do cat $i |  awk  '$3 ~ /[0-9]/'| awk '{print $1"\t"$6"\t"$7"\t"$2"\t"$3"\t"$8"\t"$4"\t"$10"\t"$5}' ; done > all_mutations
