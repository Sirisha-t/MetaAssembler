#REF=../data/strep_1Ref.fasta
#DATA=../data/simreads.fasta
#DATA=../data/test.faa
#DATA=/home/tsirisha/Projects/Project2/data/benchmark_data/cami/cami.2x.dna.fasta
DATA=/home/tsirisha/Projects/Project2/data/benchmark_data/cami/cami.2x.prot.faa
#DATA=/home/tsirisha/Projects/Project1/impp-nf/data/supp/sds1/strep.5x.fasta
#DATA=/home/tsirisha/Projects/Project2/data/benchmark_data/marine/marine.2x.prot.faa
#DATA=../data/strep_2x.faa
BIN=/home/tsirisha/Projects/Project2/MM_protein/bin
## Protein run
#/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $BIN/mm_index  --seq_type prot --k 10 --w 3 --target_file $DATA --query_file $DATA --outfile marine.2x.prot.overlaps &> marine.2x.prot.log


#NUCl run
#/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $BIN/mm_index  --seq_type nucl --k 11 --w 3 --target_file $DATA --query_file $DATA --outfile strep.5x.dna.overlaps &> strep.5x.dna.log


## Nucleotide run
#/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $BIN/mm_index  --seq_type nucl --k 15 --w 3 --target_file $DATA --query_file $DATA --outfile cami.2x.k15w3.test


## Protein run
/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $BIN/mm_index  --seq_type prot --k 11 --w 3 --target_file $DATA --reduced_alphabet 1 --query_file $DATA --outfile cami.2x.prot.k11w3.dssp5 &> cami.2x.prot.dssp5.log

#DATA=/home/tsirisha/Projects/Project2/data/benchmark_data/marine/marine.10x.prot.faa
#DATA=../data/strep_2x.faa
#BIN=/home/tsirisha/Projects/Project2/MM_protein/bin
## Protein run
#/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $BIN/mm_index  --seq_type prot --k 10 --w 3 --target_file $DATA --query_file $DATA --outfile marine.10x.prot.overlaps &> marine.10x.prot.log
