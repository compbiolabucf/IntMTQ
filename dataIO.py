import numpy as np



def save_genename(G):
	with open('RNASeq_Data/GeneName_RNASeq.txt','w') as f:
		for i in range(len(G)):
			f.write(G[i])
			if i == len(G)-1:
				continue
			else:
				f.write('\n')

def save_transcriptname(T):
	with open('RNASeq_Data/TranscriptName_RNASeq.txt','w') as f:
		for i in range(len(T)):
			for j in range(len(T[i])):
				if j == len(T[i])-1:
					f.write(T[i][j])
				else:
					f.write(T[i][j])
					f.write('\t')   # separate each transcript
			if i == len(T)-1:
				continue
			else:
				f.write('\n')       # separate each gene 

def save_transcriptlength(L):
	with open('RNASeq_Data/TranscriptLength_RNASeq.txt','w') as f:
		for i in range(len(L)):
			for j in range(len(L[i])):
				if j == len(L[i])-1:
					f.write(str(L[i][j]))
				else:
					f.write(str(L[i][j]))
					f.write('\t')   # separate each transcript
			if i == len(L)-1:
				continue
			else:
				f.write('\n')       # separate each gene 

def save_samplename(S):
	with open('RNASeq_Data/SampleName_RNASeq.txt','w') as f:
		for s in range(len(S)):
			if s == len(S)-1:
				f.write(S[s])
			else:
				f.write(S[s])
				f.write('\n')

def save_CellLineNM(CellLine):
	with open('RNASeq_Data/CellLine_RNASeq.txt','w') as f:
		for i in range(len(CellLine)):
			if i == len(CellLine)-1:
				f.write(CellLine[i])
			else:
				f.write(CellLine[i])
				f.write('\n')

def save_qtable(sample_name,Q):
	with open('RNASeq_Data/qtable/%s.txt'%sample_name,'w') as f:
		for i in range(len(Q)):
			a = Q[i]
			(p,q) = a.shape
			for m in range(p):
				for n in range(q):
					if n == q-1:
						if m == p-1:
							f.write(str(a[m][n]))
						else:
							f.write(str(a[m][n]))
							f.write(';')
					else:
						f.write(str(a[m][n]))
						f.write(',')
			if i == len(Q)-1:
				continue
			else:
				f.write('\n')   # separate each gene


class LoadRNASeqData:
	'''
	read original RNASeq Data (obtained from readcounts)
	'''

	def loadCellLineNM(self,path):
		with open(path,'r') as f:
			content = f.read()
			CellLine = []
			for row in content.split('\n'):
				CellLine.append(row)
		return CellLine

	def loadsampleNM(self,path):
		with open(path,'r') as f:
			content = f.read()
			S = []
			for row in content.split('\n'):
				S.append(row)
		return S

	def loadgenelist(self,path):
		with open(path,'r') as f:
			content = f.read()
			G = []
			for row in content.split('\n'):
				G.append(row)
		return G

	def loadisolist(self,path):
		with open(path,'r') as f:
			content = f.read()
			I = []
			for row in content.split('\n'):
				isoform_NM = []
				for item in row.split('\t'):
					isoform_NM.append(item)
				I.append(isoform_NM)
		return I

	def loadisolength(self,path):
		with open(path,'r') as f:
			content = f.read()
			L = []
			for row in content.split('\n'):
				transcript_length = []
				for item in row.split('\t'):
					transcript_length.append(int(item))
				L.append(transcript_length)
		return L

	def loadfeatByiso(self,path):
		with open(path,'r') as f:
			content = f.read()
			F = []
			for row in content.split('\n'):
				if row == '':
					F.append(np.zeros((1,1)))
				elif row == '0.0':
					F.append(np.zeros((1,1)))
				else:
					ff = []
					for rr in row.split(';'):
						each_row = []
						for cc in rr.split(','):
							each_row.append(float(cc))
						ff.append(each_row)
					ff = np.asarray(ff)
					F.append(ff)
		return F
