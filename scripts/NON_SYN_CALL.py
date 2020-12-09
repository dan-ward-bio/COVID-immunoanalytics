import pandas as pd
import sys
from multiprocessing import Pool
from datetime import datetime



def DansFunction(fastadict,SeqIndexList,ReferenceKey):

	OutputDF = pd.DataFrame()

	for seqindex in SeqIndexList:

		# Determine Refence
		ReferenceNucleotide = fastadict[ReferenceKey][seqindex]


		# Obtain list of all amino acids in relative position.
		TempAminoAcid = []
		for faseq in fastadict.values():
			TempAminoAcid.append(faseq[seqindex])
		# Remove Hyphen and X values values
		CleanTempAminoAcid = [AAT for AAT in TempAminoAcid if AAT.upper() not in ['-','X']]


		PosDict = {}
		UniqueAminoAcids = list(set(CleanTempAminoAcid))
		for A in UniqueAminoAcids:
			PosDict[A] = TempAminoAcid.count(A)


		# Determine abundance of Ref Nucleotide
		RefAbundance = PosDict[ReferenceNucleotide]

		OutputRow = {'POS':seqindex+1,
					 'REF':ReferenceNucleotide,
					 'REF_FREQ':RefAbundance,
					 'REF_REL_FREQ':round((RefAbundance/len(CleanTempAminoAcid)),2)}

		ALTString = []
		ALTFreqString = []
		ALTRelString = []
		SeqNames = []

		if len(list(PosDict.keys())) > 1:
			for KA in list(PosDict.keys()):
				if KA != ReferenceNucleotide:
					ALTString.append(KA)
					ALTFreqString.append(str(int(PosDict[KA])))
					ALTRelString.append(str(round((PosDict[KA]/len(CleanTempAminoAcid)),10)))

					MutantSeqs = []
					for FaseqTwo in list(fastadict.keys()):
						if fastadict[FaseqTwo][seqindex] == KA:
							MutantSeqs.append(FaseqTwo)

					SeqNames.append('|'.join(MutantSeqs))



		OutputRow['ALT'] = ':'.join(ALTString)
		OutputRow['ALT_FREQ'] = ':'.join(ALTFreqString)
		OutputRow['ALT_REL_FREQ'] = ':'.join(ALTRelString)
		OutputRow['ALT_Sequences'] = ':'.join(SeqNames)

		OutputDF = OutputDF.append(OutputRow,ignore_index=True)


	return OutputDF




################################################################################
################################################################################
################################################################################



InputFile = sys.argv[1]
Threads = int(sys.argv[2])


UpdatedSeqID = False

# Read in the aligned fasta sequences
fastadict = {}
Identity = 0
with open(InputFile) as file_one:
	for line in file_one:
		line = line.strip()
		if not line:
			continue


		if line.startswith(">"):
			active_sequence_name = line[1:]
			if active_sequence_name in fastadict:
				Identity +=1
				UpdatedSeqID = True
				active_sequence_name = active_sequence_name+'_{}'.format(str(Identity))
				fastadict[active_sequence_name] = ''

			else:
				fastadict[active_sequence_name] = ''
			continue

		sequence = line
		fastadict[active_sequence_name] += sequence.upper()


# Extract reference sequence
ReferenceKey = list(fastadict.keys())[0]
print('Reference Sequence: {}'.format(ReferenceKey))


# Check all Alignment Seqs are Same Length
FirstSeq = True
FastaSeqLength = 0
for fastaseq in fastadict.values():
	if FirstSeq == True:
		FastaSeqLength = len(fastaseq)
		FirstSeq = False
	if len(fastaseq) != FastaSeqLength:
		sys.exit('ERROR: Alignment file: sequence length not equal')


HDFLST = list(range(FastaSeqLength))
HomoDFInputIndexBlocks = [HDFLST[i:i + 50] for i in range(0, len(HDFLST), 50)]
MTDFOVI = list(zip([fastadict]*len(HomoDFInputIndexBlocks),HomoDFInputIndexBlocks,[ReferenceKey]*len(HomoDFInputIndexBlocks)))


with Pool(processes=Threads) as pool:
	AlignedDFMultiThreadOupt = pool.starmap(DansFunction,MTDFOVI)


MergedOutput = pd.concat(AlignedDFMultiThreadOupt)
MergedOutput = MergedOutput.sort_values(by=['POS'])

MergedOutput.to_csv('{}.AminoAcidCounts.csv'.format(InputFile),index=None)

if UpdatedSeqID == True:
	UpdatedNameFasta = open('{}_Updated_Name_Seqs.fasta'.format(InputFile),'w')
	for SeqID in list(fastadict.keys()):
		UpdatedNameFasta.write('>{}\n{}\n'.format(SeqID,fastadict[SeqID]))
	UpdatedNameFasta.close()
