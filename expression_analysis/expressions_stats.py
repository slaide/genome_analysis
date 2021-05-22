#!/usr/bin/env python3

import sys
import os

for file in os.listdir('.'):
	if not file.endswith('.csv'):
		continue

	# e.g. bin_2_expression.csv
	bin_num=int(file.split('_')[1])

	with open(file) as file:
		lines=file.readlines()
		special_features=lines[-5:]
		lines=lines[:-5]
	
		num_genes=len(lines)
		
		sizes={}
		for line in lines:
			[gene_id,count]=line.split(',')
			count=int(count[:-1])
			if count in sizes:
				sizes[count]+=1
			else:
				sizes[count]=1
	
		sizes={key:sizes[key] for key in sorted(sizes)}
		print("genes expressed in bin {:2d} : {:4d} / {:4d} = {:7.4f}%".format(bin_num,num_genes-sizes[0],num_genes,(1.0-sizes[0]/num_genes)*100))
