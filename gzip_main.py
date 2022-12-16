# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
import numpy as np
from huffmantree import HuffmanTree


class GZIPHeader:
	''' class for reading and storing GZIP header fields '''

	ID1 = ID2 = CM = FLG = XFL = OS = 0
	MTIME = []
	lenMTIME = 4
	mTime = 0

	# bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
	FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0   
	
	# FLG_FTEXT --> ignored (usually 0)
	# if FLG_FEXTRA == 1
	XLEN, extraField = [], []
	lenXLEN = 2
	
	# if FLG_FNAME == 1
	fName = ''  # ends when a byte with value 0 is read
	
	# if FLG_FCOMMENT == 1
	fComment = ''   # ends when a byte with value 0 is read
		
	# if FLG_HCRC == 1
	HCRC = []
		
	input = []
	output = []
		
	
	def read(self, f):
		''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

		# ID 1 and 2: fixed values
		self.ID1 = f.read(1)[0]  
		if self.ID1 != 0x1f: return -1 # error in the header
			
		self.ID2 = f.read(1)[0]
		if self.ID2 != 0x8b: return -1 # error in the header
		
		# CM - Compression Method: must be the value 8 for deflate
		self.CM = f.read(1)[0]
		if self.CM != 0x08: return -1 # error in the header
					
		# Flags
		self.FLG = f.read(1)[0]
		
		# MTIME
		self.MTIME = [0]*self.lenMTIME
		self.mTime = 0
		for i in range(self.lenMTIME):
			self.MTIME[i] = f.read(1)[0]
			self.mTime += self.MTIME[i] << (8 * i) 				
						
		# XFL (not processed...)
		self.XFL = f.read(1)[0]
		
		# OS (not processed...)
		self.OS = f.read(1)[0]
		
		# --- Check Flags
		self.FLG_FTEXT = self.FLG & 0x01
		self.FLG_FHCRC = (self.FLG & 0x02) >> 1
		self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
		self.FLG_FNAME = (self.FLG & 0x08) >> 3
		self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4
					
		# FLG_EXTRA
		if self.FLG_FEXTRA == 1:
			# read 2 bytes XLEN + XLEN bytes de extra field
			# 1st byte: LSB, 2nd: MSB
			self.XLEN = [0]*self.lenXLEN
			self.XLEN[0] = f.read(1)[0]
			self.XLEN[1] = f.read(1)[0]
			self.xlen = self.XLEN[1] << 8 + self.XLEN[0]
			
			# read extraField and ignore its values
			self.extraField = f.read(self.xlen)
		
		def read_str_until_0(f):
			s = ''
			while True:
				c = f.read(1)[0]
				if c == 0: 
					return s
				s += chr(c)
		
		# FLG_FNAME
		if self.FLG_FNAME == 1:
			self.fName = read_str_until_0(f)
		
		# FLG_FCOMMENT
		if self.FLG_FCOMMENT == 1:
			self.fComment = read_str_until_0(f)
		
		# FLG_FHCRC (not processed...)
		if self.FLG_FHCRC == 1:
			self.HCRC = f.read(2)
			
		return 0
			



class GZIP:
	''' class for GZIP decompressing file (if compressed with deflate) '''

	gzh = None
	gzFile = ''
	fileSize = origFileSize = -1
	numBlocks = 0
	f = None
	

	bits_buffer = 0
	available_bits = 0		

	
	def __init__(self, filename):
		self.gzFile = filename
		self.f = open(filename, 'rb')
		self.f.seek(0,2)
		self.fileSize = self.f.tell()
		self.f.seek(0)



	def decompress(self):
		''' main function for decompressing the gzip file with deflate algorithm '''
		
		numBlocks = 0

		# get original file size: size of file before compression
		origFileSize = self.getOrigFileSize()
		print(origFileSize)
		
		# read GZIP header
		error = self.getHeader()
		if error != 0:
			print('Formato invalido!')
			return
		
		# show filename read from GZIP header
		print(self.gzh.fName)
		
		
		# MAIN LOOP - decode block by block
		BFINAL = 0	
		while not BFINAL == 1:	


			
			BFINAL = self.readBits(1)
							
			BTYPE = self.readBits(2)					
			if BTYPE != 2:
				print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks+1))
				return
			
									
			#--- STUDENTS --- ADD CODE HERE

			hft = HuffmanTree()
			
			HLIT, HDIST, HCLEN = self.read_block_info()

			comp = self.code_lengths(HCLEN)

			comp_list = self.huff_converter(comp,19)

			for i in range(len(comp_list)):
				if comp_list[i] != '':
					hft.addNode(comp_list[i], i, True)

			lengths_hlit = []
			pos = 0
			while (len(lengths_hlit) < HLIT + 257):
				bit = str(self.readBits(1))
				value = hft.nextNode(bit)
				while(value < 0):
					bit = str(self.readBits(1))
					value = hft.nextNode(bit)
				if value < 16:
					lengths_hlit.append(value)
				elif value == 16:
					for i in range(3 + self.readBits(2)):
						lengths_hlit.append(lengths_hlit[-1])
				elif value == 17:
					for i in range(3 + self.readBits(3)):
						lengths_hlit.append(0)
				elif value == 18:
					for i in range(11 + self.readBits(7)):
						lengths_hlit.append(0)
				hft.resetCurNode()

			print("lengths_hlit: ", lengths_hlit)
			
			lengths_hdist = []
			while(len(lengths_hdist) < HDIST + 1):
				bit = str(self.readBits(1))
				value = hft.nextNode(bit)
				while(value < 0):
					bit = str(self.readBits(1))
					value = hft.nextNode(bit)
				if value < 16:
					lengths_hdist.append(value)
				elif value == 16:
					for i in range(3 + self.readBits(2)):
						lengths_hdist.append(lengths_hdist[-1])
				elif value == 17:
					for i in range(3 + self.readBits(3)):
						lengths_hdist.append(0)
				elif value == 18:
					for i in range(11 + self.readBits(7)):
						lengths_hdist.append(0)
				hft.resetCurNode()
			

			print("lengths_hdist: ", lengths_hdist)

			hlit_huff = self.huff_converter(lengths_hlit, 288)
			hdist_huff = self.huff_converter(lengths_hdist, 32)

			print(hlit_huff)
			print(hdist_huff)
			
			hft_lits = HuffmanTree()
			hft_dist = HuffmanTree()

			for i in range(len(hlit_huff)):
				if hlit_huff[i] != '':
					hft_lits.addNode(hlit_huff[i], i, True)
			for i in range(len(hdist_huff)):
				if hdist_huff[i] != '':
					hft_dist.addNode(hdist_huff[i], i, True)

			decompressed = []
			pos = 0
			length_base = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258]
			extra_bits = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0]
			distance_base = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577]
			extra_bits_dist = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 0]

			while(True):
				bit = str(self.readBits(1))
				value = hft_lits.nextNode(bit)
				while(value < 0):
					bit = str(self.readBits(1))
					value = hft_lits.nextNode(bit)
				if value < 256:
					print("value literal: ", value)
					decompressed.append(value)
				elif value == 256:
					break
				elif value > 256:
					length = length_base[value - 257] + self.readBits(extra_bits[value - 257])
					print("value: ", value)
					print("length: ", length)
					bit = str(self.readBits(1))
					value = hft_dist.nextNode(bit)
					while(value < 0):
						bit = str(self.readBits(1))
						value = hft_dist.nextNode(bit)
					print("value dist: ", value)
					distance = distance_base[value] + self.readBits(extra_bits_dist[value])
					print("distance: ", distance)

					if distance < len(decompressed):
						for i in range(length):
							print("decompressed len:",len(decompressed))
							print("difference:",len(decompressed) - distance)
							decompressed.append(decompressed[len(decompressed) - distance])

				hft_dist.resetCurNode()
				hft_lits.resetCurNode()

			print("decompressed: ", decompressed)
			for i in range(len(decompressed)):
				decompressed[i] = "{:08b}".format(decompressed[i])
			#print("decompressed: ", decompressed)

			output = ""
			'''
			with open("Doc4 - byteStream.txt", "r") as f:
				linhas = f.readlines()
				for i in decompressed:
					for linha in linhas:
						if linha != "\n":
							par_linha = linha.split("    ")
							byte = par_linha[0]
							simbolo = par_linha[1][1]
							if byte == i:
								output += simbolo
			'''
								
			
			print("output: ", output)
			# update number of blocks read
			numBlocks += 1
			print("numBlocks =", numBlocks)
		

		
		# close file			
		
		self.f.close()	
		print("End: %d block(s) analyzed." % numBlocks)

	#Ex 1
	def read_block_info(self):
		HLIT = self.readBits(5)
		HDIST = self.readBits(5)
		HCLEN = self.readBits(4)
		print("HCLEN:", HCLEN, "HLIT:", HLIT, "HDIST:", HDIST)
		return HLIT, HDIST, HCLEN

	#Ex 2
	def code_lengths(self, HCLEN):
		comp = np.zeros(19, dtype=int)
		for i in range(HCLEN+4):
			alfabeto_comp = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
			value = self.readBits(3)
			comp[alfabeto_comp[i]] = value

		print("comp:",comp)
		return comp

	#Ex 3
	def huff_converter(self, comp, size):
		bl_count = np.zeros(size, dtype=int)

		for i in comp:
			bl_count[int(i)] += 1

		print("bl_count:",bl_count)

		code = 0
		bl_count[0] = 0
		next_code = np.zeros(size, dtype=int)

		for bits in range(1, size):
			code = (code + bl_count[bits-1]) << 1
			next_code[bits] = code

		print("next_code:",next_code)

		comp_list = []

		for i in range(len(comp)):
			leng = comp[i]
			if leng != 0:
				comp[i] = next_code[leng]
				comp_list.append('{:03b}'.format(next_code[leng]))
				next_code[leng] += 1
			else:
				comp_list.append('')
		
		print("comp final:",comp)
		print("comp list:",comp_list)
		return comp_list

	#Ex 4
	def read_lit_len(self, HLIT):
		lit_len = np.zeros(HLIT+257, dtype=int)
		for i in range(HLIT+257):
			self.readBits(1)
		print("lit_len:",lit_len)
		return lit_len	

		
	

	def getOrigFileSize(self):
		''' reads file size of original file (before compression) - ISIZE '''
		
		# saves current position of file pointer
		fp = self.f.tell()
		
		# jumps to end-4 position
		self.f.seek(self.fileSize-4)
		
		# reads the last 4 bytes (LITTLE ENDIAN)
		sz = 0
		for i in range(4): 
			sz += self.f.read(1)[0] << (8*i)
		
		# restores file pointer to its original position
		self.f.seek(fp)
		
		return sz		
	

	def getHeader(self):  
		''' reads GZIP header'''

		self.gzh = GZIPHeader()
		header_error = self.gzh.read(self.f)
		return header_error
		

	def readBits(self, n, keep=False):
		''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

		while n > self.available_bits:
			self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
			self.available_bits += 8
		
		mask = (2**n)-1
		value = self.bits_buffer & mask

		if not keep:
			self.bits_buffer >>= n
			self.available_bits -= n

		return value
	
	def read_block_format(self):
		HLIT = self.readBits(5)
		HDIST = self.readBits(5)
		HCLEN = self.readBits(4)
		return HLIT, HDIST, HCLEN

def search_bit_by_bit(buffer, hft, verbose=False):

		lv = 0
		l = len(buffer)
		terminate = False
		code = ""

		while not terminate and lv < l:
			
			nextBit = buffer[lv]
			code = code + nextBit
			
			pos = hft.nextNode(nextBit)
						
			if pos != -2:
				terminate = True
			else:
				lv = lv + 1

		if verbose:
			if pos == -1:
				print("Code '" + buffer + "' not found!!!")
			elif pos == -2:
				print("Code '" + buffer + "': not found but prefix!!!")
			else:
				print("Code '" + buffer + "' found, alphabet position: " + str(pos) )

		return pos


if __name__ == '__main__':

	# gets filename from command line if provided
	fileName = "FAQ.txt.gz"
	if len(sys.argv) > 1:
		fileName = sys.argv[1]			

	# decompress file
	gz = GZIP(fileName)
	gz.decompress()
	