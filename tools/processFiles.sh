#!/bin/bash

# Get the file
file=$1 
# Get the filename without the extension
filename="${file%.*}" 
# Paste the '.processed.gtf' extension to the filename
processed="$filename.processed.gtf"
# Paste the '.bed' extension to the filename 
final="$filename.bed"

echo "Processing File " $file

echo "Keeping only protein coding genes"
grep -E 'protein' $file > $processed

echo "Converting to .bed file"
gtf2bed < $processed > $final
