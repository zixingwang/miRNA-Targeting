## "The characterization of microRNA-mediated gene regulation as impacted by both target site location and seed match type"
## By Zixing Wang and Wenlong Xu

import sys, string, os, math, time, re
from rnahybrid import RNAHybrid
from rnafold import RNAFold
###############################################################################
## count the seed match number of different types in different regions,
## and the total free energy of each mRNA and its potential target
class Targets_cal:
	def __init__ (self, home_dir):
		self._home=home_dir
		self._rna_hybrid = RNAHybrid(self._home+'/RNAhybrid')
		self._rna_fold = RNAFold(self._home+'/RNAfold')
		self._cut_off_delta_G = -15.0
		self._cut_off_delta_delta_G = -10.0
	def compute_delg(mrna_seq,mirna_seq):
		delg = 0.0
		rna_hybrid=self._rna_hybrid
		output = rna_hybrid(mrna_seq, mirna_seq)
		output = output.strip()
		output = output.split(":")
		if not output or not output[0]:
			break
		delg = float(output[4])
        return delg
	def compute_mfe(mrna_seq):
		mfe_open = 0.0
		rna_fold = self._rna_fold
		output = rna_fold(mrna_seq)
		output = output.split('\n')
		if len(output) < 2:
			print output
			raise Exception("Couldn't parse output")

		m = re.search(r'\(\s*([-0-9]*[.]{0,1}[0-9]+)\)$', output[1])
		if m:
			mfe_open = float(m.group(1))
		else:
			print output
			raise Exception("Couldn't parse output")
		return mfe_open
	
	def get_miRNA_sequences(self, _mirns):
		## readin the miRNAs sequnences
		miRNASeq = {}   
		fr = open(self._home + _mirns,'r')
		aLine = fr.readline()
		while aLine != '':
			if aLine[:1] == '>':
				miRNA = aLine[1:-1]
				Seq = fr.readline()[:-1]
			miRNASeq[miRNA] = Seq
			aLine = fr.readline()
		fr.close()
		return miRNASeq
    def get_weight(self, _wts):		
		W = {}
		fr = open(self._home + _wts,'r')
		aLine = fr.readline()
		for region in ['up2k','5UTR','CDS','3UTR']:
			aLine = fr.readline()
			aLine = aLine[:-1].split('\t')
			W[region] = {}
			W[region]['T1'] = float(aLine[1])
			W[region]['T2'] = float(aLine[2])
			W[region]['T3'] = float(aLine[3])
			W[region]['T4'] = float(aLine[4])
			W[region]['T5'] = float(aLine[5])
		fr.close()
		return W
		
	def get_mrna_sequence(self, fdirPaht):
			## readin the ENSP sequnences
		#print 'loading mRNAs sequences ...'
		#fdirPath = './mRNAs_samples/'
		fdirs = os.listdir(fdirPath)
		ENSPSeq = {} ##readin the ENSP sequnence (promoter/up2k, 5'UTR, CDS, 3'UTR, together in a tube)
		for fdir in fdirs:
			fr = open(fdirPath+fdir,'r')
			aLine = fr.readline()
			while aLine != '':
				aLine = aLine[:-1].split('\t')
				ENSP = aLine[0]
				Seq = aLine[1]
				if ENSP not in ENSPSeq:
					ENSPSeq[ENSP] = {'up2k':'', '5UTR':'', '3UTR':'', 'CDS':''}
				region = fdir.replace('hsa_','').replace('_ENSP.seqs','')
				ENSPSeq[ENSP][region] = Seq
				aLine = fr.readline()
			fr.close()
		return ENSPSeq
	
		

if __name__ == '__main__':
   # sc = SparkContext(appName="SparkDao Test")
	dir_curr = os.getwd()
	targets = Targets_cal(dir_curr)
	miRNAsets = targets.get_miRNA_sequence(dir_curr+'/miRNAs_samples/samples.fa')
	ENSPs = targets.get_mrna_sequences(dir_curr+'/mRNAs_samples/')
	W = targets.get_weight(dir_curr + '/weight.txt')
	if '' in ENSPs:
		ENSPs.remove('')
    
	for miRNA in miRNAsets:
		print 'computing', miRNA, '...'
		miR_seq = miRNASeq[miRNA]

		fw = open(fdirPathW + miRNA.replace('*','-star') + '.txt', 'w')
		fw.write('miRNA\tENSP\tGene')
		fw.write('\tPromoters-2t8A1\tPromoters-2t8\tPromoters-2t7A1\tPromoters-2t7\tPromoters-1t8GU')
		fw.write('\t5\'UTRs-2t8A1\t5\'UTRs-2t8\t5\'UTRs-2t7A1\t5\'UTRs-2t7\t5\'UTRs-1t8GU')
		fw.write('\tCDS-2t8A1\tCDS-2t8\tCDS-2t7A1\tCDS-2t7\tCDS-1t8GU')
		fw.write('\t3\'UTRs-2t8A1\t3\'UTRs-2t8\t3\'UTRs-2t7A1\t3\'UTRs-2t7\t3\'UTRs-1t8GU')
		fw.write('\tTotalDeltaEnergyOfDuplex\tTotalDeltaDeltaEnergy\n')

		for ENSP in ENSPs:
			print 'computing', miRNA, 'with', ENSP, '...'
			miR_TG = miRNA.replace('*','-star') +'\t' + ENSP+'\t' + SYM[ENSP]

			TotalDeltaG = 0
			TotalDeltaGs = 0
			F = {}
			
			for region in ['up2k','5UTR','CDS','3UTR']:
				miR_TG_seq = ENSPSeq[ENSP][region]
				F[region] = [0 for x in range(5)]

				if miR_TG_seq:
					seed1t8 = seq_module.compensate(seq_module.reverse(miR_seq[0:8]))
					seed2t8A1 = seed1t8[0:7] + 'A'
					seed2t8 = seed1t8[0:7]
					seed2t7A1 = seed1t8[1:7] + 'A'
					seed2t7 = seed1t8[1:7]

					TG_len = len(miR_TG_seq)
					lastp = 0
					for p in range(15,TG_len-20):
						dG2 = 0.0
						dGs = 0.0
						cons_score_seed = 0.0
						
						if re.finditer(seed2t8A1, miR_TG_seq[p:p+8]) and p-lastp > 8:
							
							dG2 = targets.compute_delg(miR_TG_seq[p-30:p+28], miR_seq)
							dGs = targets.compute_mfe( miR_TG_seq[p-30:p+28] )
							if dG2 < 0 :
								F[region][0] += 1
								TotalDeltaG += dG2 * W[region]['T1']
								TotalDeltaGs += dGs * W[region]['T1']                  
							lastp = p
							continue

						if re.finditer(seed2t8, miR_TG_seq[p:p+8]) and p-lastp > 8:
							dG2 = targets.compute_delg(miR_TG_seq[p-30:p+28], miR_seq)
							dGs = targets.compute_mfe( miR_TG_seq[p-30:p+28] )
							if dG2 < 0 :
								F[region][0] += 1
								TotalDeltaG += dG2 * W[region]['T2']
								TotalDeltaGs += dGs * W[region]['T2']                  
							lastp = p
							continue

						if re.finditer(seed2t7A1, miR_TG_seq[p:p+8]) and p-lastp > 8:
							dG2 = targets.compute_delg(miR_TG_seq[p-31:p+27], miR_seq)
							dGs = targets.compute_mfe( miR_TG_seq[p-31:p+27] )
							if dG2 < 0 :
								F[region][2] += 1
								TotalDeltaG += dG2 * W[region]['T3']
								TotalDeltaGs += dGs * W[region]['T3']                  
							lastp = p
							continue

						if re.finditer(seed2t7, miR_TG_seq[p:p+8]) and p-lastp > 8:
							dG2 = targets.compute_delg(miR_TG_seq[p-31:p+27], miR_seq)
							dGs = targets.compute_mfe( miR_TG_seq[p-31:p+27] )
							if dG2 < 0 :
								F[region][2] += 1
								TotalDeltaG += dG2 * W[region]['T4']
								TotalDeltaGs += dGs * W[region]['T4']                  
							lastp = p
							continue
						
						M = match_module.GU( seq_module.reverse(miR_seq[0:8]), miR_TG_seq[p:p+8] )
						if M[0] == 7 and M[1] == 1 and M[2] == 0 and p-lastp > 8:
							dG2 = targets.compute_delg(miR_TG_seq[p-30:p+28], miR_seq)
							dGs = targets.compute_mfe( miR_TG_seq[p-30:p+28] )
							if dG2 < 0 :
								F[region][4] += 1
								TotalDeltaG += dG2 * W[region]['T5']
								TotalDeltaGs += dGs * W[region]['T5']                  
							lastp = p
							
			TotalDeltaDeltaG = TotalDeltaG - TotalDeltaGs
			if TotalDeltaG < cut_off_delta_G and TotalDeltaDeltaG < cut_off_delta_delta_G:
				fw.write(miR_TG)
				for region in ['up2k','5UTR','CDS','3UTR']:
					for j in range(5):
						fw.write('\t'+ str(F[region][j]))
				fw.write('\t' + str(TotalDeltaG))
				fw.write('\t' + str(TotalDeltaDeltaG))
				fw.write('\n')
			
		fw.close()



