import re
import sys
import string
from Bio.Seq import Seq


## DNA codon table for bacteria and archae ( genetic code : 11)
aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
			  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
			  "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
			  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
			  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
			  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
			  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
			  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
			  "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
			  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
			  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}


def translate(offset, seq):
	out_sixfr = ""
	for x in range(int(offset), len(seq), 3):
		codon = seq[x:x + 3]
		if len(codon) < 3:
			break
		elif "N" in codon:
			out_sixfr += 'X'
		else:
			if codon in aa_dict:
				out_sixfr += aa_dict[codon]
			else:
				out_sixfr += 'X'
	return out_sixfr

## Reads input fasta reads file and converts read names to ids
def  read_fasta(fasta_file):
	read_map = {}
	seq_map = {}
	all_reads = set()
	counter = 0
	with open(fasta_file, 'r') as fasta:
		for line in fasta:
			if line[0] == '>':
				rname = line[1:line.rfind('/')]
				read_map[rname] = counter
				all_reads.add(counter)
			
			else:
				seq_map[counter] = line
				counter += 1
			

	return all_reads, read_map, seq_map


## Get reads mapping to predicted edges and paths using bwa sam file mapping info
def get_predicted_reads(bwa_name, name_map):
	predicted_set = set()
	with open(bwa_name, 'r') as bwa:
		for line in bwa:
			if line.startswith('@'):
				continue
			fields = line.strip().split('\t')	
			read_id = name_map[fields[0]]

			flag = fields[2]
			if flag != '*':
				predicted_set.add(read_id)
	
	return predicted_set

## Get reads predicted by FGS
def get_fgs_predictions(gff_name, name_map):
	predicted_set = set()
	with open(gff_name, 'r') as gff:
		for line in gff:
			if line.startswith('#'):
				continue
			fields = line.strip().split('\t')	
			rname = line[1:line.rfind('/')]
			read_id = name_map[fields[0][0:fields[0].rfind('\t')-1]]
			predicted_set.add(read_id)
	
	return predicted_set


## Six frame translate reads not predicted after all rounds of gene-calling 
def translate_unpredicted_reads(all_unpredicted_reads, name_map, seq_map, outfile):
	with open(outfile, 'w') as out:
		for read in all_unpredicted_reads:
			seq = seq_map[read].upper().strip()
			seq_rc = Seq(seq)
			seq_rc = str(seq_rc.reverse_complement().strip())

			# Create a dictionary to store the 6-frame translation
			translation = {"+1": "","+2": "","+3": "","-1": "","-2": "","-3": ""}
		
			for i in range(3):
				frame = i
				translation[f"+{frame+1}"] = translate(frame,seq)
				translation[f"-{frame+1}"] = translate(frame,seq_rc)

			longest_frame = max(translation, key=lambda k: len(translation[k]))
			
			if '*' in translation[longest_frame] and len(translation[longest_frame]) > 20:
				prefix = translation[longest_frame].split('*')[0]
				out.write(f">{read}_{longest_frame}\n{seq_map[read]}")
			elif len(translation[longest_frame]) > 20:
				out.write(f">{read}_{longest_frame}\n{seq_map[read]}")

			'''
			for frame, seq in translation.items():
				if '*' in seq:
					prefix = seq.split('*')[0]
					if len(prefix) > 20:
						out.write(f">{read}_{frame}\n{seq_map[read]}")
				else:
					if len(seq) > 20:
						out.write(f">{read}_{frame}\n{seq_map[read]}")
			'''

## Write all predicted reads from fgs, edges and paths to file
def output_predicted_reads(all_predicted_reads, name_map, seq_map, outfile):
	with open(outfile, 'w') as out:
		for read in all_predicted_reads:
			out.write(f">{read}\n")
			out.write(f"{seq_map[read]}")



if __name__ == "__main__":
	if len(sys.argv) == 6:
		read_gff_name  = sys.argv[1]
		bwa_edgename = sys.argv[2]
		bwa_pathname = sys.argv[3]
		fasta_read = sys.argv[4]
		out_dirname = sys.argv[5]
	else:
		print("ERROR: Please make sure you enter correct input arguments")
		exit(1)
		
	name_map = {}
	seq_map = {}
	edge_set = set()
	path_set = set()
	read_set = set()
	all_reads = set()

	## Read input fasta file and store ids and seqs
	all_reads, name_map, seq_map = read_fasta(fasta_read)

	## Functions to reads all the predictions at each stage -- fgs, edge and path
	edge_set = get_predicted_reads(bwa_edgename, name_map)
	path_set = get_predicted_reads(bwa_pathname, name_map)
	read_set = get_fgs_predictions(read_gff_name, name_map)


	## getting union of all predictions
	all_predicted_reads = read_set.union(edge_set, path_set)
	reads_to_be_translated = all_reads - all_predicted_reads


	seqout =  out_dirname+'/reads.translated.fasta'
	predout = out_dirname+'/reads.filtered.fasta'
	mapout = out_dirname+'/reads.namemap.txt'
	## Functions to translate the unpredicted reads in all 6 frames and output predicted and translated reads
	translate_unpredicted_reads(reads_to_be_translated, name_map, seq_map, seqout)
	output_predicted_reads(all_predicted_reads, name_map, seq_map, predout)

	#with open(mapout, 'w') as out:
	#	for key, value in name_map.items():
	#		out.write(f"{key}\t{value}\n")

