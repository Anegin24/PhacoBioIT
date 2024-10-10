#!/bin/bash
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

# Assign the input directory provided as an argument
inputdirectory="$1"

# Define other paths and directories
ref="$script_dir/Ref/hg38chr.fa"
refcnvkit="$script_dir/Ref/veriseqPGSver2.cnn"

# Create necessary directories
mkdir -p "$inputdirectory/bam"
mkdir -p "$inputdirectory/processbam"
mkdir -p "$inputdirectory/finalbam"
mkdir -p "$inputdirectory/result"

# Retrieve fastq files from the input directory
forwardread=($inputdirectory/fastq/*_R1*)

if [ ${#forwardread[@]} -eq 0 ]; then
    echo "No fastq files found in $inputdirectory. Exiting."
    exit 1
fi

# Iterate through fastq files and perform processing
for forwardread in "${forwardread[@]}"; do
reverseread=$(echo $forwardread | sed 's\R1\R2\g')
    # Extract the output filename from the fastq file name
    outputfilename=$(basename "$forwardread" | cut -f1 -d_)

bwa mem -t 24 $ref $forwardread $reverseread | samtools view -h -b -o $inputdirectory/bam/${outputfilename}.bam  

cd $inputdirectory/processbam

samtools view -b -F 0xc -@ 24 $inputdirectory/bam/${outputfilename}.bam -o ${outputfilename}.filtered.bam

samtools sort -@ 24 -n ${outputfilename}.filtered.bam -o ${outputfilename}.sorted.n.bam

samtools fixmate -@ 24 -m ${outputfilename}.sorted.n.bam ${outputfilename}.fixmate.bam

samtools sort -@ 24 ${outputfilename}.fixmate.bam -o ${outputfilename}.sorted.p.bam

samtools markdup -r -@ 24 ${outputfilename}.sorted.p.bam $inputdirectory/finalbam/${outputfilename}.dedup.bam

cd $inputdirectory/result 
mkdir -p "${outputfilename}"
cd $inputdirectory/result/${outputfilename}
echo -e 'Number-of-total-reads\tNumber-of-mapped-reads\tNumber-of-reads-after-filtering' > $inputdirectory/result/${outputfilename}/${outputfilename}.QC.tsv
echo -e "$(samtools view -c $inputdirectory/bam/${outputfilename}.bam)\t$(samtools view -c $inputdirectory/processbam/${outputfilename}.filtered.bam)\t$(samtools view -c $inputdirectory/finalbam/${outputfilename}.dedup.bam)" >> $inputdirectory/result/${outputfilename}/${outputfilename}.QC.tsv

rm $inputdirectory/bam/*.bam
rm $inputdirectory/processbam/*.bam

/home/anegin97/miniconda3/bin/cnvkit.py batch -m wgs "$inputdirectory/finalbam/${outputfilename}.dedup.bam" -r "$refcnvkit" -p 24 -d "$inputdirectory/result/${outputfilename}/"
cd "$inputdirectory/result/${outputfilename}"
zip "${outputfilename}.dedup.zip" "${outputfilename}.dedup.cnr" "${outputfilename}.dedup.cns" "${outputfilename}.dedup.call.cns" "${outputfilename}.QC.tsv" "${outputfilename}"
mv $inputdirectory/result/${outputfilename}/*.zip $inputdirectory/result
rm -r $inputdirectory/result/${outputfilename}
done

