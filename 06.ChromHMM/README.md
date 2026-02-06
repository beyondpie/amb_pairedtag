# Runnning ChromHMM
### Zhaoning Wang, test 1
### 2024-09-14

## Read this: https://www.encodeproject.org/documents/d0a10470-b049-4da1-9de2-01449ddfa6a5/@@download/attachment/ChromHMM_tutorial.pdf

## 1. Convert bam files to bed files using bedtools bamtobed

```
bedtools bamtobed -i 326_OPC_NN.H3K27ac.sorted.bam.rmpileup.sorted.bam > ./bed/326_OPC_NN.H3K27ac.bed &
bedtools bamtobed -i 326_OPC_NN.H3K27me3.sorted.bam.rmpileup.sorted.bam > ./bed/326_OPC_NN.H3K27me3.bed &
bedtools bamtobed -i 326_OPC_NN.H3K4me1.sorted.bam.rmpileup.sorted.bam > ./bed/326_OPC_NN.H3K4me1.bed &
bedtools bamtobed -i 326_OPC_NN.H3K9me3.sorted.bam.rmpileup.sorted.bam > ./bed/326_OPC_NN.H3K9me3.bed &


bedtools bamtobed -i 327_Oligo_NN.H3K27me3.sorted.bam.rmpileup.sorted.bam > ./bed/327_Oligo_NN.H3K27me3.bed &
bedtools bamtobed -i 327_Oligo_NN.H3K4me1.sorted.bam.rmpileup.sorted.bam > ./bed/327_Oligo_NN.H3K4me1.bed &
bedtools bamtobed -i 327_Oligo_NN.H3K27ac.sorted.bam.rmpileup.sorted.bam > ./bed/327_Oligo_NN.H3K27ac.bed &
bedtools bamtobed -i 327_Oligo_NN.H3K9me3.sorted.bam.rmpileup.sorted.bam > ./bed/327_Oligo_NN.H3K9me3.bed &
```



## 2. Create an empty file called cellmarkfiletable.txt in ~/data.

```
326_OPC_NN H3K27ac	326_OPC_NN.H3K27ac.bed
326_OPC_NN H3K27me3	326_OPC_NN.H3K27me3.bed
326_OPC_NN H3K4me1	326_OPC_NN.H3K4me1.bed
326_OPC_NN H3K9me3	326_OPC_NN.H3K9me3.bed
326_OPC_NN ATAC		326_OPC_NN.ATAC.bed
```

### OR:

```
326_OPC_NN      H3K27ac		326_OPC_NN.H3K27ac.bed
326_OPC_NN      H3K27me3	326_OPC_NN.H3K27me3.bed
326_OPC_NN      H3K4me1		326_OPC_NN.H3K4me1.bed
326_OPC_NN      H3K9me3		326_OPC_NN.H3K9me3.bed
327_Oligo_NN    H3K27ac		327_Oligo_NN.H3K27ac.bed
327_Oligo_NN    H3K27me3	327_Oligo_NN.H3K27me3.bed
327_Oligo_NN    H3K4me1		327_Oligo_NN.H3K4me1.bed
327_Oligo_NN    H3K9me3		327_Oligo_NN.H3K9me3.bed
```

## 3. Binarize the tracks converted to .bed with this command:

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /mnt/c/ChromHMM/CHROMSIZES/mm10.txt /mnt/d/projects/amb_PT_bam/NN_test/bed /mnt/d/projects/amb_PT_bam/NN_test/bed/cellmarkfiletable.txt /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData
```

### 5 mark version

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /mnt/c/ChromHMM/CHROMSIZES/mm10.txt /mnt/d/projects/amb_PT_bam/NN_test/bed /mnt/d/projects/amb_PT_bam/NN_test/bed/cell5markfiletable.txt /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData5marks
```

### In Chenxu's paper, binsize = 1000 bp. I should try that! Not very good to look at.

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar BinarizeBed -b 1000 /mnt/c/ChromHMM/CHROMSIZES/mm10.txt /mnt/d/projects/amb_PT_bam/NN_test/bed /mnt/d/projects/amb_PT_bam/NN_test/bed/cellmarkfiletable.txt /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData1000
```

### 2 cell types:

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /mnt/c/ChromHMM/CHROMSIZES/mm10.txt /mnt/d/projects/amb_PT_bam/NN_test/bed /mnt/d/projects/amb_PT_bam/NN_test/bed/morecelltable.txt /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData_2types
```

## 4. Learn a model:

```
java -mx64000M -jar /mnt/c/ChromHMM/ChromHMM.jar LearnModel /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output 10 mm10
```

### 5 mark version

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar LearnModel /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData5marks /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks 10 mm10
```

### 2 cell types

```
java -mx100000M -jar /mnt/c/ChromHMM/ChromHMM.jar LearnModel /mnt/d/projects/amb_PT_bam/NN_test/bed/binarizedData_2types /mnt/d/projects/amb_PT_bam/NN_test/bed/twocelltypes_output 10 mm10
```


## 6. Split bed file by state

```
#!/bin/bash

# Base name of your BED file
input_bed="326_OPC_NN_10_segments.bed"

# Loop over E1 to E10 using a traditional for loop
for i in $(seq 1 10)
do
    entry="E$i"
    # Use awk to filter based on the 4th column and output to a new file
    awk -v entry="$entry" '$4 == entry' "$input_bed" > "${entry}.bed"
    echo "${entry}.bed generated."
done
```



## 7. Plot heatmap

```
computeMatrix reference-point -S /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K27ac.sorted.bam.filtered.smooth300.bw \
							     /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K27me3.sorted.bam.filtered.smooth300.bw \
								 /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K4me1.sorted.bam.filtered.smooth300.bw \
								 /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K9me3.sorted.bam.filtered.smooth1000.bw \
							  -R /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E1.bed \
							     /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E2.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E3.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E4.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E5.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E6.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E7.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E8.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E9.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E10.bed \
							  -a 2500 \
							  -b 2500 \
							  -p max \
							  --missingDataAsZero \
							  -out /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.matrix 

plotHeatmap --matrixFile /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.matrix  \
			--outFileName /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.png \
			--colorMap Reds Oranges Greens Purples
		



computeMatrix reference-point -S /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.bw \
								 /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K27ac.sorted.bam.filtered.smooth300.bw \
							     /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K27me3.sorted.bam.filtered.smooth300.bw \
								 /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K4me1.sorted.bam.filtered.smooth300.bw \
								 /mnt/d/projects/amb_PT_bam/NN_test/326_OPC_NN.H3K9me3.sorted.bam.filtered.smooth1000.bw \
							  -R /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E1.bed \
							     /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E2.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E3.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E4.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E5.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E6.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E7.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E8.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E9.bed \
								 /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/E10.bed \
							  -a 2500 \
							  -b 2500 \
							  -p max \
							  --missingDataAsZero \
							  -out /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.withATAC.matrix 

plotHeatmap --matrixFile /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.withATAC.matrix  \
			--outFileName /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.withATAC.png \
			--colorMap Greys Reds Oranges Greens Purples
			

plotHeatmap --matrixFile /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.withATAC.matrix  \
			--outFileName /mnt/d/projects/amb_PT_bam/NN_test/bed/326_OPC_NN_output_5marks/split_state_bed/chromhmm.5marks.10states.withATAC.pdf \
			--colorMap Greys Reds Oranges Greens Purples
```
   
