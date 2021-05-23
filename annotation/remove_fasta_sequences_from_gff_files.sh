# gff files must not contain nucleotide sequences (fasta format/fasta section) for htseq expression analysis
for bin_name in $(ls prokka_out); do

echo $bin_name ; date
gff_file="prokka_out/$bin_name/$bin_name.gff"

# get number of lines in file
nextcmd="wc -l $gff_file"
line_count=$($nextcmd)
line_count=$(echo $line_count | awk -v delim=" " '{print substr($0, 0, index($0, delim)-1)}')
#echo $line_count

# get line number where fasta section starts
nextcmd="awk -v line='##FASTA' '$0==line{print NR}' ${gff_file}"
#fasta_header_line=$($nextcmd)
fasta_header_line=$(awk --assign=line='##FASTA' '$0==line{print NR}' $gff_file)

# count number of lines in fasta section (excluding header), hopefully
#nextcmd="sed -n '${fasta_header_line},${line_count}p' ${gff_file} | sed /\>/g | wc -l"
#echo $nextcmd
fasta_section_line_count=$(sed -n ${fasta_header_line},${line_count}p ${gff_file} | sed /\>/g | wc -l)
#echo $fasta_section_line_count

# + 1 to include the header line
if ! test $fasta_section_line_count -eq $(($line_count - $fasta_header_line + 1))
then
echo "data after fasta section. that sucks"
fi

last_gff_line=$(($fasta_header_line-1))
nextcmd="$(head -${last_gff_line} ${gff_file} > prokka_out/${bin_name}/${bin_name}_nofa.gff)"
echo $nextcmd
$nextcmd

done
