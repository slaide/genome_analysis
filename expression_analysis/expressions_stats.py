#!/usr/bin/env python3

import sys
import os
import re
import matplotlib.pyplot as plt

bin_stats={}

for expression_level_file in os.listdir('.'):
	if not expression_level_file.endswith('.csv'):
		continue

	bin_num=int(expression_level_file.split('_')[1])

	annotation_filename="../annotation/prokka_out/bin_{}/bin_{}_nofa.gff".format(bin_num,bin_num)
	with open(annotation_filename) as file:
		lines=file.readlines()
		#remove lines beginning with hash and remove end of line characters
		lines=[re.sub('[\r\n]','',line) for line in lines if line[0]!='#']
		#then remove empty lines
		lines=[line for line in lines if len(line)>0]

		segments=[line.split('\t') for line in lines]
		segments=[(sequence_name,source,feature,start,end,noclue1,noclue2,noclue3,features) for (sequence_name,source,feature,start,end,noclue1,noclue2,noclue3,features) in segments if feature=="CDS"]

		sequence_data={}

		for (sequence_name,source,feature,start,end,noclue1,noclue2,noclue3,attributes) in segments:
			attributes=attributes.split(';')
			named_attributes={}

			for feature in attributes:
				(name,value)=feature.split('=')
				named_attributes[name]=value

			#ignore subsequenct occurences of annotated regions within the same sequence
			gene_id=named_attributes["ID"]
			if gene_id not in sequence_data:
				data={}
				data["product"]=named_attributes["product"]
				#copy these attributes, if available
				for possible_attribute in "gene Name".split():
					if possible_attribute in named_attributes:
						data[possible_attribute]=named_attributes[possible_attribute]

				sequence_data[gene_id]=data
			else:
				assert False, "this is bad"

		bin_stats[bin_num]=sequence_data
		

	with open(expression_level_file) as file:
		lines=file.readlines()
		special_lines=lines[-5:]
		lines=lines[:-5]
	
		num_genes=len(lines)
		
		sizes={}
		for line in lines:
			[gene_id,count]=line.split(',')
			count=int(count[:-1])

			sequence_data[gene_id]['count']=count

	stats=bin_stats[bin_num]

	hypothetical_proteins_expressed=0
	hypothetical_proteins_not_expressed=0
	for gene_id in stats.keys():
		gene=stats[gene_id]
		product_value=gene["product"]
		if product_value=="hypothetical protein":
			if gene['count']==0:
				hypothetical_proteins_not_expressed+=1
			else:
				hypothetical_proteins_expressed+=1

for bin_key in bin_stats.keys():
	expression_levels=[
		gene['count']
		for gene
		in bin_stats[bin_key].values()
	]
	expression_level_stuff={}
	for level in expression_levels:
		if level in expression_level_stuff:
			expression_level_stuff[level]+=1
		else:
			expression_level_stuff[level]=1

	expression_level_stuff=[els for els in expression_level_stuff.items()]
	#sum of all levels to normalize to max expression per bin
	whole_expression=sum([level for (level,value) in expression_level_stuff])
	#sum of all genes to normalize to number of genes per bin
	whole_count=sum([value for (level,value) in expression_level_stuff])

	#allow all genes, not only those expressed (expression level >0) and those with non-unique expression level (count>1). the latter is used to make the plot more compact, and the former is used to restrict the plot to 'really' expressed genes
	allow_all=False
	#level=intensity of expression
	#count=number of genes at given expression level
	xv=[level for (level,count) in expression_level_stuff if allow_all or (level>0 and count >1)]
	yv=[count for (level,count) in expression_level_stuff if allow_all or (level>0 and count >1)]

	plt.scatter(xv,yv,label="{}".format(bin_key))
plt.title("Genes with same expression level per bin\nexcluding genes with unique expression level\nor expression level of 0")
plt.ylabel("number of genes")
plt.xlabel("expression level")
plt.legend()
plt.show()

# calculate fraction of hypothetical proteins
fraction_hypothetical_proteins=[
	(
		sum([
			1
			for gene
			in filter(lambda gene: gene['product']=="hypothetical protein",bin_stats[bin_key].values())
		])/len(bin_stats[bin_key].keys()),
		bin_key
	)
	for bin_key
	in bin_stats.keys()
]
print(max([value for value,num in fraction_hypothetical_proteins]))
print(min([value for value,num in fraction_hypothetical_proteins]))
print(sum([value for value,num in fraction_hypothetical_proteins])/len(fraction_hypothetical_proteins))