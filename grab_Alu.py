#!/usr/bin/env python

with open("Alu_sublibrary","U") as input_file:
	tokens = input_file.readlines()
	start = False
	cur_seq = ""
	cur_alu = ""
	with open("Alu_seqs.tsv","w") as output_file:
		for index,line in enumerate(tokens):
			# Found an Alu element to start
			if line.startswith(";ORIGIN"):
				start = True
				cur_alu = "> " + tokens[index+1] # Add next line as index
				cur_alu_mark = tokens[index+1]
			elif start == False:
				pass
			elif line.startswith(cur_alu_mark):
				print "skipping header"
				pass
			elif line.startswith(";LOCUS"):
				print "writing to file"
				output_file.write(cur_alu + cur_seq[:-1] + "\n")
				start = False
				cur_seq = ""
			elif start == True:
				# Add sequence
				cur_seq += line.strip()
