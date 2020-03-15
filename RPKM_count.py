def fileread1(filename):
	"""
		input:bedtools coverage result file
		chrom   start   end     peakname overlap_readcounts	length1	length2	precentage
		chr1    4784062 4784416 peak_1 32      354     354     1.0000000 
		output:dict[peakname]=overlap_readcounts
		outdict[peak_1]=32
	"""
	f=open(filename,'r')
	peaknamelist = []
	valuelist = {}
	for str_x in f:
		list_x = str_x[:-1].split('\t')
		valuelist[list_x[3]]=list_x[4]
	return valuelist
def fileread2(filename):
	"""
		input:peaks.bed
		chrom	start	end	peakname
		chr1    4784062 4784416 peak_1
		output:
		outlist = [[chr1,4784062,4784416,peak_1]]
	"""
	f = open(filename,'r')
	outlist = []
	for str_x in f:
		list_x = str_x[:-1].split('\t')
		outlist.append(list_x)
	return outlist

def resultwrite(peaknamelist,valuelist,allread):
	"""
		if overlap_readcounts==0 or readlength ===0:
			RPKM = 0
		else:
			RPKM = (overlap_readcounts*1000000*1000)/(total_readcounts*readlength)
	"""
	for peakname in peaknamelist:
		try:
			readcount = float(valuelist[peakname[3]])
		except:
			readcount = 0
		readlength = float(peakname[2])-float(peakname[1])
		if readlength == 0:
			RPKM=0
		else:
			RPKM = (readcount*1000000*1000)/(allread*readlength)
		print('\t'.join(['\t'.join(peakname),str(RPKM)]))
if __name__=='__main__':
	import sys
	filename1 = sys.argv[1]              # bedtools coverage result file
	filename2 = sys.argv[3]              # peaks.bed file
	allread = float(sys.argv[2])         # total read of ChIP-seq sample
	L1=fileread2(filename2)
	L2=fileread1(filename1)
	resultwrite(L1,L2,allread) 
