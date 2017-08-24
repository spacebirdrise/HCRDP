import pybedtools
import csv
import sys
import os 
class DataParser:
	''' This is the main class of the tool.  From here, you can call various methods that work to convert and analyze genomic data files:
		User can choose to dump all info into this class call and call each function in this order: exonToIntron(), intronExtender(), then MainParser() or call each individually with different arguments
		Calling each function in order will yield the files desired if the string arguments in the constructor are all defined.
		Whenever you call a function make sure that the following file names are empty or do not exist:
		When calling exonToIntron: 1 file is saved: intron_file.bed
		When calling intronExtender: 1 file is saved: extended_intron_file.bed
		When calling mainParser: 39 files are saved: extended_lifted_mouse_circRNA_file.bed hcf_elmcf.bed, hcf_elmcf_same_start.bed, hcf_elmcf_same_end.bed, hcf_elmcfss_sine.bed, hcf_elmcfse_sine.bed, hcf_elmcfsss_unextended.bed, hcf_elmcfses_unextended.bed, hesu_nodups.bed, heeu_nodups.bed, hcb_sine.bed, hcbs_nodups.bed, hcbs_reextended.bed, hc_extended.bed, mc_same.bed, mc_same_sine.bed, mcss_nodups.bed, forced_liftover_mcss.bed, flm_start_extended.bed, flm_end_extended.bed, fse_b1b2.bed, fee_b1b2.bed, fseb_unextended.bed, feeb_unextended.bed, mcb_both.bed, mcbb_nodups.bed, introns_mcbb.bed, imcbb_unextended.bed, forced_liftover_mcf_human.bed, humanCircRNAfinalextended.bed, hcf_normal.bed, hcfn_nodups.bed, hcrpm.bed, hcrpm_nodups.bed, cofmv.bed, comhvp.bed, comparison_of_mouse_human_final.bed, cofmvv_use.bed, narrow_list_human_mouse.bed, nlhm_final.bed   
		(it might be a good idea to set up a separate empty directory prior to caling these methods to contain these files)
		The parameters hcf, mcf, mclf, hrsinef, and mrb1b2f, must be defined to use this code.
		If you want to look at another genome, you must call the class again to redefine elements from the new genome.
		The other parameters here can be defined later or redefined in calls to the functions, and the function definition of the parameters take priority.
		This is made to work with Python v2.7.12.
		
		:param ef: the name of the file containing the exons of the genome of interest in bed format 	
		:type ef: string
		:param hcf: the name of the file containing human circular RNA data in bed format
		:type hcf: string
		:param mcf: the name of the file containing the circular RNA data for the genome of interest in bed format
		:type mcf: string
		:param mclf: the name of the file containing the circular RNA data for the genome of interest lifted over to human coordinates using the liftOver tool from the genome browser
		:type mclf: string
		:param hrsinef: the name of the file containing the data for the SINEs from human repeats in the genome in bed format
		:type hrsinef: string 
		:param mrb1b2f: the name of the file containing the SINE equivalents contained the repeats of the genome of interest
		:type mrb1b2f: string
		:param extend_sine: the number of nt that is the max distance away from the circRNA for a SINE to be considered flanking to a circRNA (default 2000)
		:type extend_sine: int
		:param extend_circRNA: the max number of nucleotides that can be considered a negligible distance when determining whether the circRNA in the genome of interest is considered to be in the same position as in the human genome (default 50)
		:type extend_circRNA: int
		:param extend_intron: the amount of nucleotides to extend introns to act as a buffer to overlap circRNAs in order to be considered flanking to the circRNA (default 10)
		:type extend_intron: int
		:param comp_distance_buffer_high: the higher limit of how many nucleotides apart the start end coordinates of a circRNA in a given genome must be from the start and end coordinates of the other genome in order to be considered corresponding circRNAs, used in creating the nlhm_final file (default 50)
		:type comp_distance_buffer_high: int
		:param comp_distance_buffer_low: the lower limit for the number of nuleotides apart the start and end coordinates of a circRNA in a given genome must be from the start and end coordinates of the other genome in order to be considered corresponding circRNAs, used in creating the nlhm_final file (default -50)
		:type comp_distance_buffer_low: int
		'''
	def __init__(self, ef, hcf, mcf, mclf, hrsinef, mrb1b2f, extend_sine, extend_circRNA, extend_intron, comp_distance_buffer_high, comp_distance_buffer_low):
		''' This is the constructor for the DataParser class '''
		self.ef = ef
		try:
			self.hcf = hcf
			self.mcf = mcf
			self.mclf = mclf
			self.hrsinef = hrsinef	
			self.mrb1b2f = mrb1b2f
			if hcf == null:
				raise ValueError('hcf was set to null')
			

			if mcf == null:
				raise ValueError('mcf was set to null')
			

			if mclf == null:
				raise ValueError('mclf was set to null')
				
				
			if hrsinef == null:
				raise ValueError('hrsinef was set to null')
			

			if mrb1b2f == null:
				raise ValueError('mrb1b2f was set to null')
				
					
		except ValueError:
			print "Must define hcf, mcf, mclf, hrsinef, and mrb1b2f for code to run."
			
			
			
		
		if self.extend_sine == null:
			self.extend_sine = 2000
			
			
			
		self.extend_circRNA = extend_circRNA
		if self.extend_circRNA == null:
			self.extend_circRNA = 50
			
			
			
		self.extend_intron = extend_intron
		if self.extend_intron == null:
			self.extend_intron = 10
			
			
		self.comp_distance_buffer_high = comp_distance_buffer_high
		if self.comp_distance_buffer_high == null:
			self.comp_distance_buffer_high = 50
			
			
		self.comp_distance_buffer_low = comp_distance_buffer_low
		if self.comp_distance_buffer_low == null:
			self.comp_distance_buffer_low = -50
			
			
			
	def exonToIntron(self, ef2):
		'''Converts an exon file to an intron file 
		
		:param ef2: string of the bed file containing information on the exons of the genome of interest (default self.ef)
		:type ef: string
		'''
		if ef2 != null:
			ef = ef2
			
			
		if ef == null:
			print "Must define exon file in either call to function or initial class call"
			return
			
			
		c1 = list()
		with open(ef) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					c1.append(word)
		names = c1[0::12]
		chroms = c1[1::12]
		starts = c1[8::12]
		ends = c1[9::12]
		starts2 = list()
	#have to remove the first element of every starting exon list because the beginning start exon for a given gene is the edge of the gene
		for x in starts:
			one_start = x.split(",")
			one_start.remove(one_start[0])
			starts2.append(one_start)

		starts3 = list()
	#have to remove any excess spaces that are at the end any of the elements, so they do not get extracted with the introns coordinates
		for x in starts2:
			one_start = x
			if one_start != []:
				if one_start[-1] == '':
					one_start.remove(one_start[-1])
				starts3.append(one_start)

		starts4 = list()
	#dump all the individual coords into a large list
		for x in starts3:
			for y in x:    
				starts4.append(y)

		ends2 = list()
	#remove all of the extra empty spaces at the end and the last number of each list
		for x in ends:
			one_end = x.split(",")
			if(one_end[-1] == ''):
				one_end.remove(one_end[-1])
			one_end.remove(one_end[-1])
			if one_end != []:
				ends2.append(one_end)

		ends3 = list() 
	#dump all of the coords int one list
		for x in ends2:
			for y in x:
				ends3.append(y)

		f = open('intron_file.bed', 'w')		
		k = 0
	#ends come before starts because the ends of the exon are the beginnings of the introns and the starts of the exons are the ends of the introns
		for i in range(len(names)):
			for j in range(len(starts3[i])):
				print >> f, chroms[i] + "\t" + ends3[k] + "\t" + starts4[k] + "\t" + names[i]
				k = k + 1

		f.close()
		print "The file containing the introns based on the exon file you submitted are saved under intron_file.bed in the current directory"
	def intronExtender(self, inf2):
		''' This extends the intron coordinates by 10 nt in both directions 
		
			:param inf2: Optional string parameter that represents an intron file that the user may want to extend separately (default "intron_file.bed", the file produced at the end of the last exonToIntron() call
			:type inf: string
		'''
		
		if inf2 != null:
			inf = inf2
			
			
		if inf == null:
			inf = "intron_file.bed"
			
		if inf == null:
			print "Must define intron_file.bed in either call to function or initial call to class, or run a successful call to exonToIntron()"
			return
			
			
		c1 = list()
		with open(inf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					c1.append(word)
		chroms = c1[0::4]
		starts = c1[1::4]
		ends = c1[2::4]
		names = c1[3::4]
		starts2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start - extend_intron 
			starts2.append(one_start)
		ends2 = list()
		for x in ends:
			one_end = int(x)
			one_end = one_end + extend_intron
			ends2.append(one_end)
		f = open('extended_intron_file.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i]

		f.close()
	def mainParser(self, feiom, fextend_sine, fextend_circRNA, fextend_intron, fcomp_distance_buffer_high, fcomp_distance_buffer_low):
		'''The main function with the purpose of analyzing, comparing, and producing files with respect to the genome of interest's relation to the human genome.
		
		:param feiom: a string representing the file containing the extended introns of the genome of interest(default eiom as definedby the class) (default2 "extended_intron_file.bed")
		:type feiom: string
		:param fextend_sine: the number of nt that is the max distance away from the circRNA for a SINE to be considered flanking to a circRNA (default 2000) 
		:type fextend_sine: int
		:param fextend_circRNA: the max number of nucleotides that can be considered a negligible distance when determining whether the circRNA in the genome of interest is considered to be in the same position as in the human genome (default extend_circRNA) 
		:type fextend_circRNA: int
		:param fextend_intron: the amount of nucleotides to extend introns to act as a buffer to overlap circRNAs in order to be considered flanking to the circRNA (default 10)
		:type fextend_intron: int
		:param fcomp_distance_buffer_high: the higher limit of how many nucleotides apart the start end coordinates of a circRNA in a given genome must be from the start and end coordinates of the other genome in order to be considered corresponding circRNAs, used in creating the nlhm_final file (default 50)
		:type fcomp_distance_buffer_high: int
		:param fcomp_distance_buffer_low: the lower limit for the number of nuleotides apart the start and end coordinates of a circRNA in a given genome must be from the start and end coordinates of the other genome in order to be considered corresponding circRNAs, used in creating the nlhm_final file (default -50)
		:type fcomp_distance_buffer_low: int
	
		'''
		if feiom != null:
			eiom = feiom
		
		
		if eiom == null:
			eiom = "extended_intron_file.bed"
			
			
		if eiom == null:
			print "Must define extended_intron_file.bed in either class or function call, or run a successful intronExtender"
			return

		if fextend_sine != null:
			extend_sine = fextend_sine
		
		
		if fextend_circRNA != null:
			extend_circRNA = fextend_circRNA
			
			
		if fextend_intron != null:
			extend_intron = fextend_intron
			
			
		if fcomp_distance_buffer_high != null:
			comp_distance_buffer_high = fcomp_distance_buffer_high
			
			
		if fcomp_distance_buffer_low != null:
			comp_distance_buffer_low = fcomp_distance_buffer_low
			
		database = list()
		with open(mclf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		ends2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start - extend_amount
			starts2.append(one_start)


		for x in ends:
			one_end = int(x)
			one_end = one_end + extend_amount
			ends2.append(one_end)


		f = open('extended_lifted_mouse_circRNA_file.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()

		#now it's ready to intersect
		hc = pybedtools.BedTool(hcf)
		mcl = pybedtools.BedTool('extended_lifted_mouse_circRNA_file.bed')
		f = open('hcf_elmcf.bed', 'w')
		print >> f, hc.intersect(mcl, wa=True)
		f.close()
		#have to extend each side one at a time
		#code to extend hc_samev2.bed by 2000 on start side
		database = list()
		with open('hcf_elmcf.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start - extend_sine
			starts2.append(one_start)


		f = open('hcf_elmcf_same_start.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + ends[i] + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#Now do the end side
		database = list()
		with open('hcf_elmcf.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		ends2 = list()
		for x in ends:
			one_end = int(x)
			one_end = one_end + extend_sine
			ends2.append(one_end)


		f = open('hcf_elmcf_same_end.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + starts[i] + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#intersect to get files for both sides
		hcse = pybedtools.BedTool('hcf_elmcf_same_start.bed')
		hrsine = pybedtools.BedTool(hrsinef)
		f = open('hcf_elmcfss_sine.bed', 'w')
		print >> f, hcse.intersect(hrsine, wa=True)
		f.close()
		#intersect #2
		hcse = pybedtools.BedTool('hcf_elmcf_same_end.bed')
		hrsine = pybedtools.BedTool(hrsinef)
		f = open('hcf_elmcfse_sine.bed', 'w')
		print >> f, hcse.intersect(hrsine, wa=True)
		f.close()
		#unextend each file 
		database = list()
		with open('hcf_elmcfss_sine.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start + extend_sine
			starts2.append(one_start)


		f = open('hcf_elmcfsss_unextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + ends[i] + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#unextend the ends
		database = list()
		with open('hcf_elmcfse_sine.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		ends2 = list()
		for x in ends:
			one_end = int(x)
			one_end = one_end - extend_sine
			ends2.append(one_end)

		f = open('hcf_elmcfses_unextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + starts[i] + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#duplicates
		os.system("cat 'hcf_elmcfsss_unextended.bed' | sort | uniq > hesu_nodups.bed")
		os.system("cat 'hcf_elmcfses_unextended.bed' | sort | uniq > heeu_nodups.bed")
		#now have to intersect the two files to get a listof introns that have both requirements
		hs = pybedtools.BedTool('hesu_nodups.bed')
		he = pybedtools.BedTool('heeu_nodups.bed')
		f = open('hcb_sine.bed', 'w')
		print >> f, hs.intersect(he, wa=True)
		f.close()
		#delete duplicates
		os.system("cat 'hcb_sine.bed' | sort | uniq > hcbs_nodups.bed")
		#must be reextended by 50 on both sides
		database = list()
		with open('hcbs_nodups.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		ends2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start - extend_circRNA
			starts2.append(one_start)


		for x in ends:
			one_end = int(x)
			one_end = one_end + extend_circRNA
			ends2.append(one_end)


		f = open('hcbs_reextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#have to remake mc_samev3
		database = list()
		with open(hcf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		ends2 = list()
		for x in starts:
			one_start = x
			if one_start !='chromStart':
				one_start = int(x)
				one_start = one_start - extend_circRNA
			starts2.append(one_start)


		for x in ends:
			one_end = x
			if one_end != 'chromEnd':
				one_end = int(x)
				one_end = one_end + extend_circRNA
			ends2.append(one_end)



		f = open('hc_extended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		mcsl = pybedtools.BedTool(mclf)
		hcsr = pybedtools.BedTool('hc_extended.bed')
		f = open('mc_same.bed', 'w')
		print >> f, mcsl.intersect(hcsr, wa=True)
		f.close()
		#now mc_samev3(lifted) can be intersected with hcs_reextended with mcsl.intersect(hcsr, wa=True)
		mcsl = pybedtools.BedTool('mc_same.bed')
		hcsr = pybedtools.BedTool('hcbs_reextended.bed')
		f = open('mc_same_sine.bed', 'w')
		print >> f, mcsl.intersect(hcsr, wa=True)
		f.close()
		#duplicates remove
		os.system("cat 'mc_same_sine.bed' | sort | uniq > mcss_nodups.bed")
		#force liftover of mc_same
		data = list()
		with open('mcss_nodups.bed') as tsv:
			 for line in csv.reader(tsv, dialect="excel-tab"):
					 for word in line:
							 data.append(word)

		 
		chroms = data[0::12]
		starts = data[1::12]
		ends = data[2::12]
		names = data[3::12]
		data1 = data[4::12]
		data2 = data[5::12]
		data3 = data[6::12]
		data4 = data[7::12]
		data5 = data[8::12]
		data6 = data[9::12]
		data7 = data[10::12]
		data8 = data[11::12]
		database = list()
		with open(mcf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms2 = database[0::12]
		starts2 = database[1::12]
		ends2 = database[2::12]
		names2 = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]

		f = open('forced_liftover_mcss.bed', 'w')

		for i in range(len(chroms)):
			for j in range(len(chroms2)):
				if names2[j] == names[i]:
					print >> f, chroms2[j] + "\t" + starts2[j] + "\t" + ends2[j] + "\t" + names2[j] + "\t" + datab1[j] + "\t" + datab2[j] + "\t" + datab3[j] + "\t" +  datab4[j] + "\t" + datab5[j] + "\t" + datab6[j] + "\t" + datab7[j] + "\t" + datab8[j]


		f.close()

		#forced_liftover_mcss.bed now is all of the mousecircRNA with the same coords as humancircRNA that have sines within 2000 nt that are in the coords of mouse genome
		#all that's left now is to use the mouseB1 and B2 file to determine which of these circRNAs have B1 or B2 within 2000 nt
		# start with an extension of 2000 nt in both directions of the current mcss bed file
		database = list()
		with open('forced_liftover_mcss.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()

		for x in starts:
			one_start = int(x)
			one_start = one_start - extend_sine
			starts2.append(one_start)

		f = open('flm_start_extended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + ends[i] + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]

		f.close()

		#ENDS
		database = list()
		with open('forced_liftover_mcss.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		ends2 = list()
		for x in ends:
			one_end = int(x)
			one_end = one_end + extend_sine
			ends2.append(one_end)

		f = open('flm_end_extended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + starts[i] + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#intersect each with B1B2
		fse = pybedtools.BedTool('flm_start_extended.bed')
		mrb1b2 = pybedtools.BedTool(mrb1b2f)
		f = open('fse_b1b2.bed', 'w')
		print >> f, fse.intersect(mrb1b2, wa=True)
		f.close()
		#ends
		fee = pybedtools.BedTool('flm_end_extended.bed')
		mrb1b2 = pybedtools.BedTool(mrb1b2f)
		f = open('fee_b1b2.bed', 'w')
		print >> f, fee.intersect(mrb1b2, wa=True)
		f.close()
		#unextend the two
		database = list()
		with open('fse_b1b2.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()

		for x in starts:
			one_start = int(x)
			one_start = one_start + extend_sine
			starts2.append(one_start)

		f = open('fseb_unextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + ends[i] + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#ends
		database = list()
		with open('fee_b1b2.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)



		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		ends2 = list()

		for x in ends:
			one_end = int(x)
			one_end = one_end - extend_sine
			ends2.append(one_end)

		f = open('feeb_unextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + starts[i] + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#intersect the two
		import pybedtools
		fse = pybedtools.BedTool('fseb_unextended.bed')
		fee = pybedtools.BedTool('feeb_unextended.bed')
		f = open('mcb_both.bed', 'w')
		print >> f, fse.intersect(fee, wa=True)
		f.close()
		#more duplicates
		os.system("cat 'mcb_both.bed' | sort | uniq > mcbb_nodups.bed")
		#intersect the bed file with the introns
		mcbb = pybedtools.BedTool('mcbb_nodups.bed')
		eiom4 = pybedtools.BedTool(eiom)
		f = open('introns_mcbb.bed', 'w')
		print >> f, eiom4.intersect(mcbb, wa=True)
		f.close()
		#unextend the intron files
		database = list()
		with open('introns_mcbb.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms = database[0::4]
		starts = database[1::4]
		ends = database[2::4]
		names = database[3::4]
		starts2 = list()
		ends2 = list()
		for x in starts:
			one_start = int(x)
			one_start = one_start + extend_intron
			starts2.append(one_start)


		for x in ends:
			one_end = int(x)
			one_end = one_end - extend_intron
			ends2.append(one_end)


		f = open('imcbb_unextended.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i]


		f.close()
		#toget human circRNA from mouse circRNA
		data = list()
		with open('mcbb_nodups.bed') as tsv:
			 for line in csv.reader(tsv, dialect="excel-tab"):
					 for word in line:
							 data.append(word)

		 
		chroms = data[0::12]
		starts = data[1::12]
		ends = data[2::12]
		names = data[3::12]
		data1 = data[4::12]
		data2 = data[5::12]
		data3 = data[6::12]
		data4 = data[7::12]
		data5 = data[8::12]
		data6 = data[9::12]
		data7 = data[10::12]
		data8 = data[11::12]
		database = list()
		with open(mclf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms2 = database[0::12]
		starts2 = database[1::12]
		ends2 = database[2::12]
		names2 = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]

		f = open('forced_liftover_mcf_human.bed', 'w')

		for i in range(len(chroms)):
			for j in range(len(chroms2)):
				if names2[j] == names[i]:
					print >> f, chroms2[j] + "\t" + starts2[j] + "\t" + ends2[j] + "\t" + names2[j] + "\t" + datab1[j] + "\t" + datab2[j] + "\t" + datab3[j] + "\t" +  datab4[j] + "\t" + datab5[j] + "\t" + datab6[j] + "\t" + datab7[j] + "\t" + datab8[j]


		f.close()


		#intersect them with human circRNA(extended since that is the consideration)
		import pybedtools
		hce = pybedtools.BedTool('hc_extended.bed')
		mcfl = pybedtools.BedTool('forced_liftover_mcf_human.bed')
		f = open('humanCircRNAfinalextended.bed', 'w')
		print >> f, hce.intersect(mcfl, wa=True, wb=True)
		f.close()

		#unextend human circRNAs
		database = list()
		with open('humanCircRNAfinalextended.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		chroms = database[0::12]
		starts = database[1::12]
		ends = database[2::12]
		names = database[3::12]
		datab1 = database[4::12]
		datab2 = database[5::12]
		datab3 = database[6::12]
		datab4 = database[7::12]
		datab5 = database[8::12]
		datab6 = database[9::12]
		datab7 = database[10::12]
		datab8 = database[11::12]
		starts2 = list()
		ends2 = list()
		for x in starts:
			one_start = x
			if one_start !='chromStart':
				one_start = int(x)
				one_start = one_start + extend_circRNA
			starts2.append(one_start)


		for x in ends:
			one_end = x
			if one_end != 'chromEnd':
				one_end = int(x)
				one_end = one_end - extend_circRNA
			ends2.append(one_end)


		f = open('hcf_normal.bed', 'w')
		for i in range(len(chroms)):
			print >> f, chroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + names[i] + "\t" + datab1[i] + "\t" + datab2[i] + "\t" + datab3[i] + "\t" + datab4[i] + "\t" + datab5[i] + "\t" + datab6[i] + "\t" + datab7[i] + "\t" + datab8[i]


		f.close()
		#more duplicates
		os.system("cat 'hcf_normal.bed' | sort | uniq > hcfn_nodups.bed")

		# get the data lined up side by side
		hce = pybedtools.BedTool('hc_extended.bed')
		mcfl = pybedtools.BedTool('forced_liftover_mcf_human.bed')
		f = open('hcrpm.bed', 'w')
		print >> f, mcfl.intersect(hce, wa=True, wb=True)
		f.close()
		#duplicates
		os.system("cat 'hcrpm.bed' | sort | uniq > hcrpm_nodups.bed")
		#code to make the file/ breakpoint
		database = list()
		with open('hcrpm_nodups.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::24]
		mstarts = database[1::24]
		mends = database[2::24]
		mnames = database[3::24]
		hchroms = database[12::24]
		hstarts = database[13::24]
		hends = database[14::24]
		hnames = database[15::24]
		f = open('cofmv.bed', 'w')
		for i in range(len(mchroms)):
			print >> f, mchroms[i] + "\t" + mstarts[i] + "\t" + mends[i] + "\t" + mnames[i] + "\t" + hchroms[i] + "\t" + hstarts[i] + "\t" + hends[i] + "\t" + hnames[i]


		f.close()
		#have to liftOver mouse back to mouse genome
		database = list()
		with open('cofmv.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::8]
		mstarts = database[1::8]
		mends = database[2::8]
		mnames = database[3::8]
		hchroms = database[4::8]
		hstarts = database[5::8]
		hends = database[6::8]
		hnames = database[7::8]
		data = list()
		with open(mcf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					data.append(word)


		chroms2 = data[0::12]
		starts2 = data[1::12]
		ends2 = data[2::12]
		names2 = data[3::12]
		f = open('comhvp.bed', 'w')
		for i in range(len(mchroms)):
			for j in range(len(chroms2)):
				if names2[j] == mnames[i]:
					print >> f, chroms2[j] + "\t" + starts2[j] + "\t" + ends2[j] + "\t" + names2[j] + "\t" + hchroms[i] + "\t" + hstarts[i] + "\t" + hends[i] + "\t" +  hnames[i]


		f.close()
		#have to unextend the RNAs on the human side
		database = list()
		with open('comhvp.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::8]
		mstarts = database[1::8]
		mends = database[2::8]
		mnames = database[3::8]
		hchroms = database[4::8]
		hstarts = database[5::8]
		hends = database[6::8]
		hnames = database[7::8]
		starts2 = list()
		ends2 = list()
		for x in hstarts:
			one_start = x
			if one_start !='chromStart':
				one_start = int(x)
				one_start = one_start + extend_circRNA
			starts2.append(one_start)


		for x in hends:
			one_end = x
			if one_end != 'chromEnd':
				one_end = int(x)
				one_end = one_end - extend_circRNA
			ends2.append(one_end)



		f = open('comparison_of_mouse_human_final.bed', 'w')
		for i in range(len(mchroms)):
			print >> f, mchroms[i] + "\t" + mstarts[i] + "\t" + mends[i] + "\t" + mnames[i] + "\t" + hchroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + hnames[i]


		f.close()
		#have to unextend human coords first
		database = list()
		with open('cofmv.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::8]
		mstarts = database[1::8]
		mends = database[2::8]
		mnames = database[3::8]
		hchroms = database[4::8]
		hstarts = database[5::8]
		hends = database[6::8]
		hnames = database[7::8]
		starts2 = list()
		ends2 = list()
		for x in hstarts:
			one_start = x
			if one_start !='chromStart':
				one_start = int(x)
				one_start = one_start + extend_circRNA
			starts2.append(one_start)


		for x in hends:
			one_end = x
			if one_end != 'chromEnd':
				one_end = int(x)
				one_end = one_end - extend_circRNA
			ends2.append(one_end)



		f = open('cofmvv_use.bed', 'w')
		for i in range(len(mchroms)):
			print >> f, mchroms[i] + "\t" + mstarts[i] + "\t" + mends[i] + "\t" + mnames[i] + "\t" + hchroms[i] + "\t" + repr(starts2[i]) + "\t" + repr(ends2[i]) + "\t" + hnames[i]


		f.close()
		#have to make sure human coords in comparison file are within 50 on both start and end of mouse
		database = list()
		with open('cofmvv_use.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::8]
		mstarts = database[1::8]
		mends = database[2::8]
		mnames = database[3::8]
		hchroms = database[4::8]
		hstarts = database[5::8]
		hends = database[6::8]
		hnames = database[7::8]
		f = open('narrow_list_human_mouse.bed', 'w')
		for i in range(len(mchroms)):
			ms = int(mstarts[i])
			me = int(mends[i])
			hs = int(hstarts[i])
			he = int(hends[i])
			x = ms - hs
			y = me - he
			if x >= distance_buffer_low:
				if x <= distance_buffer_high:
					if y >= distance_buffer_low:
						if y <= distance_buffer_high:
							print >> f, mchroms[i] + "\t" + mstarts[i] + "\t" + mends[i] + "\t" + mnames[i] + "\t" + hchroms[i] + "\t" + hstarts[i] + "\t" + hends[i] + "\t" + hnames[i]


		f.close()
		#force liftover mouse
		database = list()
		with open('narrow_list_human_mouse.bed') as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					database.append(word)


		mchroms = database[0::8]
		mstarts = database[1::8]
		mends = database[2::8]
		mnames = database[3::8]
		hchroms = database[4::8]
		hstarts = database[5::8]
		hends = database[6::8]
		hnames = database[7::8]
		data = list()
		with open(mcf) as tsv:
			for line in csv.reader(tsv, dialect="excel-tab"):
				for word in line:
					data.append(word)


		chroms2 = data[0::12]
		starts2 = data[1::12]
		ends2 = data[2::12]
		names2 = data[3::12]
		f = open('nlhm_final.bed', 'w')
		for i in range(len(mchroms)):
			for j in range(len(chroms2)):
				if names2[j] == mnames[i]:
					print >> f, chroms2[j] + "\t" + starts2[j] + "\t" + ends2[j] + "\t" + names2[j] + "\t" + hchroms[i] + "\t" + hstarts[i] + "\t" + hends[i] + "\t" +  hnames[i]


		f.close()
		print "This function has saved 39 files to your computer, but four of them are of interest:"
		print "mcbb_nodups.bed contains the circRNA from the genome of interest that correspond to circRNAs in human and contain sine equivalents on both sides within the specified sine buffer while also containing sines on both sides within the sine buffer on its human equivalent"
		print "imcbb_unextended.bed contains the flanking introns of the circRNA contained in the mcbb_nodups.bed file"
		print "hcfn_nodups.bed is the human circRNA that corresponds to the mcbb_nodups.bed circRNA"
		print "nlhm_final.bed is the bed file containing both the human circRNA and the circRNA in mcbb_nodups.bed that corresponds side by side in a bed formatted list"