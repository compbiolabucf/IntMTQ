# ========================================================
#
# count reads aligned to each gene from aligned bam file
#
# ========================================================


import os
import numpy as np 
import pandas as pd 
from dataIO import save_genename, save_transcriptname, save_transcriptlength, save_samplename, save_CellLineNM, save_qtable



def annotation(path):
	'''
	This function is to load hg19 refseq annotation with filtering the following information.
	1) isoforms in the same gene but in different chromsone or strand
	2) isoforms in chr_xxxxxx (uncertain position)
	3) non-coding isoforms (NR_xxxxxx)
	4) duplicate isoforms with the same transcript names only keep one

	Args:
	path: directory of hg19 annotation file (hg19_2018June18)

	Return:
	GeneNames,TranscriptNames,Chr,TranscriptStart,TranscriptEnd,CodingStart,CodingEnd,ExonStart,ExonEnd,ExonNum,Strand
	* 50470 transcripts after filtering, 19278 genes in total

	'''
	GeneNames = []
	TranscriptNames = []
	Chr = []
	TranscriptStart = []
	TranscriptEnd = []
	CodingStart = []
	CodingEnd = []
	ExonStart = []
	ExonEnd = []
	ExonNum = []
	Strand = []

	with open(path,'r') as f:
		content = f.read()
		for row in content.split('\n')[1:]:
			if row != '':
				refseq = row.split('\t')
				TranscriptNames.append(refseq[1])
				Chr.append(refseq[2])
				Strand.append(refseq[3])
				TranscriptStart.append(int(refseq[4]))
				TranscriptEnd.append(int(refseq[5]))
				CodingStart.append(int(refseq[6]))
				CodingEnd.append(int(refseq[7]))
				ExonNum.append(int(refseq[8]))
				ExonStart.append(refseq[9])
				ExonEnd.append(refseq[10])
				GeneNames.append(refseq[12])
	# filter isoforms in chr_xxxxxx and non-coding isoforms
	chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	for i in range(len(TranscriptNames)-1,-1,-1):
		if Chr[i] in set(chr_list) and TranscriptNames[i].split('_')[0] == 'NM':
			continue
		else:
			del TranscriptNames[i],Chr[i],Strand[i],TranscriptStart[i],TranscriptEnd[i],CodingStart[i],CodingEnd[i],ExonNum[i],ExonStart[i],ExonEnd[i],GeneNames[i]
	# filter isoforms in the same gene but in different chromsone or strand
	G = list(set(GeneNames))
	for g in G:
		IDX = [index for index,item in enumerate(GeneNames) if item == g]
		chr_g = [Chr[idx] for idx in IDX]
		strand_g = [Strand[idx] for idx in IDX]
		if len(set(chr_g))==1 and len(set(strand_g))==1:
			continue
		else:
			IDX.sort(reverse=True)
			for i in IDX:
				del TranscriptNames[i],Chr[i],Strand[i],TranscriptStart[i],TranscriptEnd[i],CodingStart[i],CodingEnd[i],ExonNum[i],ExonStart[i],ExonEnd[i],GeneNames[i]
	# filter duplicate isoforms
	duplicate_index = []
	unique_T = []
	for index,item in enumerate(TranscriptNames):
		if item not in unique_T:
			unique_T.append(item)
		else:
			duplicate_index.append(index)
	duplicate_index.sort(reverse=True)
	for i in duplicate_index:
		del TranscriptNames[i],Chr[i],Strand[i],TranscriptStart[i],TranscriptEnd[i],CodingStart[i],CodingEnd[i],ExonNum[i],ExonStart[i],ExonEnd[i],GeneNames[i]

	return GeneNames,TranscriptNames,Chr,TranscriptStart,TranscriptEnd,CodingStart,CodingEnd,ExonStart,ExonEnd,ExonNum,Strand


def build_index(path):
	'''
	Build index for aligned bam files with samtools-0.1.8.

	Args:
	path: directory of aligned bam files

	Return:
	files: list of sample names for all cancer cell lines
	CellLine: list of cancer cell line names

	'''
	files = os.listdir(path)
	CellLine = []
	for f in files:
		CellLine.append(f.split('.')[1])
		os.system('./samtools-0.1.8/samtools index bamfile/%s'%f)

	return files,CellLine


def gene_info(G,GeneNames,TranscriptNames,TranscriptStart,TranscriptEnd,Chr):
	'''
	Check start and end position of each gene. 
	Prepare isoform name and isoform length corresponding to each gene in G.

	Args:
	G: unique gene list; 
	GeneNames,TranscriptNames,TranscriptStart,TranscriptEnd,Chr (from hg19 annotation)

	Return:  (all corresponding to each gene in G)
	Gene_S: gene start position; Gene_E: gene end position; Gene_Chr: chromsone of each gene
	T: list of transcript names; L: transcript length  

	'''
	Gene_S = []
	Gene_E = []
	Gene_Chr = []
	T = []   # isoform name of each gene
	L = []   # isoform length of each gene
	for g in G:
		IDX = [index for index,item in enumerate(GeneNames) if item == g]
		a = [TranscriptStart[idx] for idx in IDX]
		b = [TranscriptEnd[idx] for idx in IDX]
		c = [Chr[idx] for idx in IDX]
		Gene_S.append(min(a))
		Gene_E.append(max(b))
		Gene_Chr.append(c[0])
		t = [TranscriptNames[idx] for idx in IDX]
		T.append(t)
		l = [abs(TranscriptEnd[idx]-TranscriptStart[idx]) for idx in IDX]
		L.append(l)

	return Gene_S,Gene_E,Gene_Chr,T,L


def parse_samfile(transcript_name):
	'''
	Parse aligned sam files for transcript 'transcript_name' (transcript_name_1.sam & transcript_name_2.sam).

	Args:
	transcript_name: current transcript name

	Return:
	intersect: pair-end reads name aligned to current transcript

	'''
	with open('bamfile/%s_1.sam'%transcript_name,'r') as f:
		sam1 = f.read()
		reads_name1 = []
		for row in sam1.split('\n'):
			if row != '':
				reads_name1.append(row)
	with open('bamfile/%s_2.sam'%transcript_name,'r') as b:
		sam2 = b.read()
		reads_name2 = []
		for row in sam2.split('\n'):
			if row != '':
				reads_name2.append(row)
	intersect = list(set(reads_name1) & set(reads_name2))
	
	return intersect


def bam2qtable(path_bamfile,path_annotation):
	'''
	Generate qtable (read counts) from aligned bam file

	Args:
	path_bamfile: directory of aligned bam files
	path_annotation: directory of annotation file (hg19_2018June18)

	Return:
	G: list of gene names
	T: list of transcript names corresponding to each gene in G
	L: list of transcript length corresponding to each gene in G
	qtable: aligned read counts, reads by isoform
	CellLine: list of cell line names

	'''

	print('build index start...')      
	files,CellLine = build_index(path_bamfile)
	print('build index end') 
	# hg19 annotation
	GeneNames,TranscriptNames,Chr,TranscriptStart,TranscriptEnd,CodingStart,CodingEnd,ExonStart,ExonEnd,ExonNum,Strand = annotation(path_annotation)
	G = list(set(GeneNames))
	Gene_S,Gene_E,Gene_Chr,T,L = gene_info(G,GeneNames,TranscriptNames,TranscriptStart,TranscriptEnd,Chr)

	save_samplename(files)
	save_CellLineNM(CellLine)
	save_genename(G)
	save_transcriptname(T)
	save_transcriptlength(L)

	for i in range(len(files)):
		print('start counting sample: %s'%files[i])     
		qtable = []      # qtable of each sample (each cancer cell line)
		for j in range(len(G)):
			reads_name = []
			for k in range(len(T[j])):
				cmd1 = './samtools-0.1.8/samtools view -f 0X0040 bamfile/%s'%files[i]+' '+'%s'%T[j][k]+' | cut -f 1'+' > '+'bamfile/%s_1.sam'%T[j][k]
				os.system(cmd1)
				cmd2 = './samtools-0.1.8/samtools view -f 0X0080 bamfile/%s'%files[i]+' '+'%s'%T[j][k]+' | cut -f 1'+' > '+'bamfile/%s_2.sam'%T[j][k]
				os.system(cmd2)
				if os.path.getsize('bamfile/%s_1.sam'%T[j][k]) == 0 or os.path.getsize('bamfile/%s_2.sam'%T[j][k]) == 0:
					reads_name.append([])
				else:
					intersect = parse_samfile(T[j][k])
					reads_name.append(intersect)
				os.system('rm bamfile/*.sam')
			union = []
			for r in reads_name:
				for rr in r:
					union.append(rr)
			if len(union) == 0:
				qtable.append(np.zeros((1,1)))
			else:
				readByiso = pd.DataFrame(data=np.zeros((len(set(union)),len(T[j]))),columns=T[j],index=list(set(union)))
				for t in range(len(T[j])):
					if len(reads_name[t]) != 0:
						readByiso.loc[reads_name[t],T[j][t]] += 1
				tmp = readByiso.copy()
				tmp[tmp>0] = 1
				new_index = []
				tmp_unique_row = []
				for index,row in tmp.iterrows():
					row_data = row.tolist()
					if row_data not in tmp_unique_row:
						tmp_unique_row.append(row_data)
						new_index.append(tmp_unique_row.index(row_data))
					else:
						new_index.append(tmp_unique_row.index(row_data))
				readByiso.index = new_index
				featByiso = pd.DataFrame(columns=T[j])
				for ix in range(len(tmp_unique_row)):
					featByiso.loc[ix] = readByiso.loc[ix].sum()
				qtable.append(featByiso.values)

		save_qtable(files[i],qtable)

	return G,T,L,qtable,CellLine


if __name__ == '__main__':

	path_bamfile = 'bamfile'
	path_annotation = 'hg19_2018June18'

	G,T,L,qtable,CellLine = bam2qtable(path_bamfile,path_annotation)