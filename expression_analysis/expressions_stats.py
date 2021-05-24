#!/usr/bin/env python3

import sys
import os
import re
import matplotlib.pyplot as plt

bin_stats={}

features_found={}

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
		
		for (sequence_name,source,feature,start,end,noclue1,noclue2,noclue3,features) in segments:
			if feature in features_found:
				features_found[feature]+=1
			else:
				features_found[feature]=1

		segments=filter(lambda line:line[2]=='CDS',segments)
		segments=list(segments)
		#print("number of segments in current bin after filtering of non-cds features:",len(segments))

		def get_product_name_from_line(line):
			attributes=line[8].split(';')
			found=None
			for attribute in attributes:
				(key,value)=attribute.split('=')
				if key=='product':
					found=value
			assert found!=None

			return found

		segments=filter(lambda line:get_product_name_from_line(line) not in "hypothetical protein,putative protein".split(','),segments)
		segments=list(segments)
		#print("number of segments in current bin after filtering of hypothetical proteins:",len(segments))

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

			if gene_id in sequence_data:
				sequence_data[gene_id]['count']=count
			else:
				# gene encodes hypothetical protein, which we ignore
				pass

	stats=bin_stats[bin_num]

	# this is now impossible, because further up hypothetical and putative proteins are removed
	hypothetical_proteins_present=False
	if hypothetical_proteins_present:
		hypothetical_proteins_expressed=0
		hypothetical_proteins_not_expressed=0
		for gene_id in stats.keys():
			gene=stats[gene_id]
			product_value=gene["product"]
			if product_value in "hypothetical protein,putative protein".split(','):
				if gene['count']==0:
					hypothetical_proteins_not_expressed+=1
				else:
					hypothetical_proteins_expressed+=1

for (key,value) in features_found.items():
	print('|',key,'|',value,'|')

# sort dictionary by bin key for better tables further down
bin_stats={key: value for key, value in sorted(bin_stats.items(), key=lambda item: item[0])}

plot_expression_levels=True

for bin_key in bin_stats.keys():
	expression_levels=filter(lambda gene:gene['product'] not in "hypothetical protein,putative protein".split(','),bin_stats[bin_key].values())
	expression_levels=sorted(expression_levels, key=lambda item: item['count'],reverse=True)
	top_n=5
	for gene in expression_levels[:top_n]:
		expression_level=gene['count']
		product_name=gene['product']
		print('|',bin_key,'|',expression_level,'|',product_name,'|')
	expression_levels=map(lambda gene:gene['count'],expression_levels)

	expression_level_abundances={}
	for level in expression_levels:
		if level in expression_level_abundances:
			expression_level_abundances[level]+=1
		else:
			expression_level_abundances[level]=1

	expression_level_abundances=list(expression_level_abundances.items())
	#sum of all levels to normalize to max expression per bin
	whole_expression=sum([level for (level,value) in expression_level_abundances])
	#sum of all genes to normalize to number of genes per bin
	whole_count=sum([value for (level,value) in expression_level_abundances])

	#allow all genes, not only those expressed (expression level >0) and those with non-unique expression level (count>1). the latter is used to make the plot more compact, and the former is used to restrict the plot to 'really' expressed genes
	allow_all=False
	#level=intensity of expression
	#count=number of genes at given expression level
	xv=[level for (level,count) in expression_level_abundances if allow_all or (level>0 and count >1)]
	yv=[count for (level,count) in expression_level_abundances if allow_all or (level>0 and count >1)]

	if plot_expression_levels:
		plt.scatter(xv,yv,label="{}".format(bin_key))

if plot_expression_levels:
	plt.title("Genes with same expression level per bin\nexcluding genes with unique expression level\nor expression level of 0")
	plt.ylabel("number of genes")
	plt.xlabel("expression level (mRNA reads mapped to gene)")
	plt.legend()
	plt.show()

# calculate fraction of hypothetical proteins
if hypothetical_proteins_present:
	fraction_hypothetical_proteins=[
		(
			sum([
				1
				for gene
				in filter(lambda gene: gene['product'] in "hypothetical protein,putative protein".split(','),bin_stats[bin_key].values())
			])/len(bin_stats[bin_key].keys()),
			bin_key
		)
		for bin_key
		in bin_stats.keys()
	]

	print_hypothetical_protein_fractions=False

	if print_hypothetical_protein_fractions:
		print(max([value for value,num in fraction_hypothetical_proteins]))
		print(min([value for value,num in fraction_hypothetical_proteins]))
		print(sum([value for value,num in fraction_hypothetical_proteins])/len(fraction_hypothetical_proteins))

# save bins that express same products
bins_per_gene={}

for bin_key in bin_stats.keys():
	# count number of genes expressing identical products
	product_abundances={}

	for gene in bin_stats[bin_key].values():
		gene_product=gene['product']

		# we have no useful information on these currently, so they are omitted for this analysis
		if gene_product not in "hypothetical protein,putative protein".split(','):
			if gene_product in bins_per_gene:
				bins_per_gene[gene_product][bin_key]=True
			else:
				bins_per_gene[gene_product]={bin_key:True}

			if gene_product not in product_abundances:
				product_abundances[gene_product]=1
			else:
				product_abundances[gene_product]+=1

	# sort dictionary for better tables further below
	product_abundances={key: value for key, value in sorted(product_abundances.items(), key=lambda item: item[1],reverse=True)}

	# print genes with copy number above threshold, or instead the n most often copied genes across the genome
	print_abundances_above_threshold_instead_of_top=True

	if  print_abundances_above_threshold_instead_of_top:
		for (product_name,product_abundance) in product_abundances.items():
			if product_abundance>5:
				print('|',bin_key,'|',product_abundance,'|',product_name,'|')
	else:
		top_n=5
		for (product_name,product_abundance) in list(product_abundances.items())[:top_n]:
			print('|',bin_key,'|',product_abundance,'|',product_name,'|')

	# count number of genes with same number of copies
	abundances={}
	for abundance in product_abundances.values():
		if abundance not in abundances:
			abundances[abundance]=1
		else:
			abundances[abundance]+=1

	# print(bin_key,abundances)

	plt.scatter([abundance for (abundance,count) in abundances.items()],[count for (abundance,count) in abundances.items()],label="{}".format(bin_key))

print("number of bins in which genes with identical product is present")
for (product,bins) in bins_per_gene.items():
	if len(bins)>13:
		print('|',len(bins), '|',product,'|')

plt.title("Number of genes with same copy number in each bin")
plt.ylabel("Number of genes with same copy number")
plt.xlabel("Copy number")
plt.legend()
# plt.show()