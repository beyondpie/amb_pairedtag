#!/bin/bash

# Define the main directory where all subclass folders are located
main_dir="/mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files"

# Loop through each subclass folder
for subclass in "$main_dir"/*; do
    # Check if it's a directory
    if [ -d "$subclass" ]; then
        echo "Processing subclass folder: $subclass"
        
        # Loop through each E*.bed file in the subclass folder
        for bed_file in "$subclass"/E*.bed; do
            # Extract the filename without the path
            filename=$(basename "$bed_file")
            
            # Subset the BED file to 10,000 peaks and overwrite the original
            shuf -n 10000 "$bed_file" > "$subclass/subset_$filename"
            
            echo "Subset created for $bed_file"
        done
    fi
done



# Directory containing the subclasses (folders)
dir="/mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files"

# Loop through each folder in the directory
for subclass in "$dir"/*/; do
  # Remove the path and trailing slash to get just the folder name
  subclass=$(basename "$subclass")
  
  # Create the necessary output directory for the subclass
  mkdir -p /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass

  # Run computeMatrix and plotHeatmap for the second set of files (4marks)
  computeMatrix reference-point -S /mnt/f/projects/amb_PT_bam/bigwigs/bigwig_ATAC_bamcoverage/$subclass.ATAC.e100.bs100.sm300.bw \
								   /mnt/f/projects/amb_PT_bam/bigwigs/bigwig_PT_DNA/$subclass.H3K27ac.e100.bs100.sm300.bw \
								   /mnt/f/projects/amb_PT_bam/bigwigs/bigwig_PT_DNA/$subclass.H3K27me3.e100.bs100.sm300.bw \
								   /mnt/f/projects/amb_PT_bam/bigwigs/bigwig_PT_DNA/$subclass.H3K4me1.e100.bs100.sm300.bw \
								   /mnt/f/projects/amb_PT_bam/bigwigs/bigwig_PT_DNA/$subclass.H3K9me3.e100.bs100.sm1000.bw \
								-R /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E1.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E2.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E3.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E4.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E5.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E6.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E7.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E8.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E9.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E10.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E11.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E12.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E13.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E14.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E15.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E16.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E17.bed \
								   /mnt/f/projects/amb_PT_bam/ChromHMM/model_bypeak_b200_s18/bed_files/$subclass/subset_E18.bed \
								-a 2500 -b 2500 -p max/2 --referencePoint center --missingDataAsZero \
								-out /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass/subset.chromhmm.5marks.18states.matrix

  # Plot heatmaps for 5 marks
  plotHeatmap --matrixFile /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass/subset.chromhmm.5marks.18states.matrix \
    --outFileName /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass/subset.chromhmm.5marks.18states.png \
    --colorMap Greys Reds Oranges Greens Purples

  plotHeatmap --matrixFile /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass/subset.chromhmm.5marks.18states.matrix \
    --outFileName /mnt/f/projects/amb_PT_bam/ChromHMM/heatmap_18states/$subclass/subset.chromhmm.5marks.18states.pdf \
    --colorMap Greys Reds Oranges Greens Purples

done

echo "Processing completed for all subclasses."


