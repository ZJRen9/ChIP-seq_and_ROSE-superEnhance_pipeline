def ROSE_annotate_fileread(filename):
	"""
		input:
		Enhancers_name	chrom	start	end	H3K27ac_number	length	IP	Input	 rank	OVERLAP_GENES   PROXIMAL_GENES  CLOSEST_GENE	PROXIMAL_gene_num	OVERLAP_gene_num
		23_Mst_H3K27Ac_peak_36002_lociStitched chr6    125404542   125505812   23  55280   281966.061  28902.458   3   Cd9	Plekhg6,Vwf Cd9 3   1
		output:
		[[chr6,125404542,125505812,23_Mst_H3K27Ac_peak_36002_lociStitched,3,"Cd9,Vwf,Plekhg6","Cd9"]]
	"""
	f =open(filename,'r')
	resultlist = []
	for str_x in f:
		if str_x[0]=="#":
			continue
		list_x = str_x[:-1].split("\t")
		peakname = list_x[0]
		chrom = list_x[1]
		start = list_x[2]
		end = list_x[3]
		rank = list_x[8]
		overlap_gene = list_x[9].split(",")
		proximal_gene = list_x[10].split(",")
		closet_gene = list_x[11].split(",")
		outlist = overlap_gene + proximal_gene + closet_gene
		outlist2 = outlist_del(outlist)
		totalgenename = ",".join(outlist2)
		if len(overlap_gene)==1:
			if overlap_gene[0]=="":
				overlap_genename="empty"
			else:
				overlap_genename=overlap_gene[0]
		else:
			overlap_genename=','.join(overlap_gene)
		resultlist.append([chrom,start,end,peakname,rank,totalgenename,overlap_genename])
		#print("\t".join([chrom,start,end,peakname,rank,totalgenename,overlap_genename]))
	return resultlist

def outlist_del(outlist):
	outlist2 = []
	for gene in outlist:
		if gene == "":
			pass
		else:
			if gene in outlist2:
				pass
			else:
				outlist2.append(gene)
	return outlist2

def RNAseq_file_to_dict(filename):
	"""
		input:
		gene_id    WT_Rep1 WT_Rep2
		0610007P14Rik  46.7721 42.5738
		0610009B22Rik  17.7982 16.178
		outdict:
		outdict[0610007P14Rik]=44.67295
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		if str_x[0]=="#":
			continue
		list_x = str_x[:-1].split("\t")
		valuelist = [float(x) for x in list_x[1:]]
		mean_value = sum(valuelist)/len(valuelist)
		genename = list_x[0]
		outdict[genename] = mean_value
	return outdict

def suitable_gene_result_write(ROSE_annotate_filelist,expdict):
	"""
		basic principle:if overlap_maxgene_exp > 10 or total_max_gene_exp < 5:suitable_gene == overlap_maxgene
		input:
		[[chr6,125404542,125505812,23_Mst_H3K27Ac_peak_36002_lociStitched,3,"Cd9,Vwf,Plekhg6","Cd9"]]
		output:
		chr6	125404542	125505812	23_Mst_H3K27Ac_peak_36002_lociStitched	3	Cd9
	"""
	for list_x in ROSE_annotate_filelist:
		genelist = list_x[5].split(",")
		overlap_genelist = list_x[6].split(",")
		maxvalue,maxgene = max_gene(genelist,expdict)
		overlap_maxvalue,overlap_maxgene = max_gene(overlap_genelist,expdict)
		most_suitable_gene = maxgene
		if (overlap_maxvalue > 10 or maxvalue < 5) and overlap_maxgene != "empty" :
			most_suitable_gene = overlap_maxgene
		print("\t".join(list_x[:5]+[most_suitable_gene]))

def max_gene(genelist,expdict):
	outlist = []
	for gene in genelist:
		try:
			v = expdict[gene]
		except:
			v = 0
		outlist.append(v)
	maxvalue = max(outlist)
	maxindex = outlist.index(maxvalue)
	maxgene = genelist[maxindex]
	#print(maxvalue,maxgene)
	return maxvalue,maxgene

if __name__=="__main__":
	import sys
	ROSE_annotate_filename = sys.argv[1]
	RNAseq_exp_filename = sys.argv[2]
	ROSE_annotate_file_list = ROSE_annotate_fileread(ROSE_annotate_filename)
	RNAseq_exp_dcit = RNAseq_file_to_dict(RNAseq_exp_filename)
	suitable_gene_result_write(ROSE_annotate_file_list,RNAseq_exp_dcit)
	
