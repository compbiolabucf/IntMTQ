# ===================================================================
# 
# This file implement the followings:
# 	1) Preprocess RNASeq Data: load RNASeq data 
# 	2) Penalized_EM function: implement method of IntMTQ with cvxpy
# 	3) Experiment E1: IntMTQ model, use RNASeq and NanoString Data, 
#	   46 cell lines & 99 genes & 304 isoforms 
#
# ===================================================================


from dataIO import LoadRNASeqData
import os, sys 
import numpy as np
import pandas as pd
import cvxpy as cvx
import math
import itertools


 
def Penalized_EM(q,s,L,maxiter,E,lambd):     
	'''
	Args:
	q: q matrix, read by isoform
	s: read counts for each read 
	L: transcript length (effective)
	E: expression of the second type of data
	lambd: hyperparameter, need to be tuned first

	Return:
	P: transcript proportion
	Expression: transcript expression

	'''

	m,n = q.shape
	P = np.ones(n)/n
	q = q/np.tile(L,(m,1))*max(L)
	M = np.diagflat([1/l for l in L])

	for i in range(maxiter):
		P_old = P
		alpha = np.sum((P_old/L)*np.sum(s))/np.sum(E)
		if math.isinf(alpha) == True:
			alpha = 0

		#### E step ####
		a = np.multiply(np.multiply(np.tile(P,(m,1)),q)/(np.tile(np.dot(q,P.T).reshape(m,1),(1,n))),np.tile(s.reshape(m,1),(1,n)))
		#### M step ####
		x = cvx.Variable(n)
		loss = -cvx.sum(a*cvx.log(x))
		reg = cvx.norm(M*x*cvx.sum(s)-alpha*E)
		objective = cvx.Minimize(loss + lambd*reg)
		constraints = [x >= 0,cvx.sum(x) == 1]
		prob = cvx.Problem(objective,constraints)
		prob.solve(solver=cvx.SCS)
		P = x.value

		P = P.T
		if np.absolute(P - P_old).max() < 1e-6:
			break
		else:
			continue
	Expression = P/L*np.sum(s)
	
	return P,Expression


def E1(G_NanoString,CellLine_NanoString,TransExpre_NanoString):
	'''
	This function implement penalized EM model, IntMTQ, with RNASeq data and NanoString data. 
	Only one penalization term has been introduced. It penalized all isoforms of each gene (99 gene).
	* RNASeq Data: 46 CellLines, 19278 Genes, 50470 Transcripts 
	* NanoString Data: 46 CellLines, 99 Genes, 304 Transcripts
	
	Args:
	NanoString data: 
		G_NanoString: list of gene names
		CellLine_NanoString: list of cell lines
		TransExpre_NanoString: NanoString expression data, transcript by CellLine

	Return:
	data_frame: transcript expression (transcript by CellLine)

	'''

	# RNASeq Data: 46 CellLine, 19278 genes, 50470 transcripts
	path_RNASeqData = 'RNASeq_Data/'
	rnaseq = LoadRNASeqData()
	SampleName_RNASeq = rnaseq.loadsampleNM(path_RNASeqData+'SampleName_RNASeq.txt')
	CellLine_RNASeq = rnaseq.loadCellLineNM(path_RNASeqData+'CellLine_RNASeq.txt')
	GeneName_RNASeq = rnaseq.loadgenelist(path_RNASeqData+'GeneName_RNASeq.txt')
	TranscriptName_RNASeq = rnaseq.loadisolist(path_RNASeqData+'TranscriptName_RNASeq.txt')
	TranscriptLength_RNASeq = rnaseq.loadisolength(path_RNASeqData+'TranscriptLength_RNASeq.txt')

	# take intersection of CellLine, Gene between RNA-Seq and NanoString 
	CellLine = list(set(CellLine_RNASeq)&set(CellLine_NanoString))  
	G = list(set(GeneName_RNASeq)&set(G_NanoString))                 
	IX1 = [GeneName_RNASeq.index(G[i]) for i in range(len(G))]

	# initialize transcript expression(Data), transcript list(Transcript)
	Data = [[] for i in range(len(CellLine))]   # CellLine by transcript
	Transcript = []

	for i in range(len(CellLine)):
		ix1 = CellLine_RNASeq.index(CellLine[i])
		qtable = rnaseq.loadfeatByiso(path_RNASeqData+'qtable/%s.txt'%SampleName_RNASeq[ix1])
		for j in range(len(G)):
			ixx1 = IX1[j]
			if i == 0:
				Transcript.append(TranscriptName_RNASeq[ixx1])
			if qtable[ixx1].shape == (1,1) and qtable[ixx1][0][0] == 0:
				for k in range(len(TranscriptName_RNASeq[ixx1])):
					Data[i].append(0)
			else:
				E = TransExpre_NanoString.loc[TranscriptName_RNASeq[ixx1],CellLine[i]].values
				s = qtable[ixx1].max(axis=1)
				q = qtable[ixx1]
				q[q>0] = 1
				L = np.asarray(TranscriptLength_RNASeq[ixx1])
				maxiter = 1000
				lambd = 10000
				P, Expression = Penalized_EM(q,s,L,maxiter,E,lambd)
				for e in Expression:
					Data[i].append(e)

	Data = np.asarray(Data)
	Transcript_flattened = list(itertools.chain(*Transcript))
	data_frame = pd.DataFrame(data=Data.T,index=Transcript_flattened,columns=CellLine)   # transcript by CellLine

	return data_frame

			


if __name__ == '__main__':

	# NanoString data
	TransExpre_NanoString = pd.read_excel(sys.argv[1])
	CellLine_NanoString = TransExpre_NanoString.columns.tolist()
	with open(sys.argv[2],'r') as f:
		G_NanoString = f.read().split('\n')

	# IntMTQ 
	data_frame = E1(G_NanoString,CellLine_NanoString,TransExpre_NanoString)
	data_frame.to_excel('E1_expression.xlsx')
