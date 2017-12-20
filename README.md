## CNV_Radar
Copy number variant (CNV) caller for detecting amplifications, deletions and copy-neutral loss of heterozygousity events in genomic data.

## Installation
This package was tested using R version 3.3.2 and bedtools v2.24.0 on the linux command line

## Dependencies 
bedtools  
R  

## Additional R Package Requirements
getopt  
data.table  

## Generate a normal cohort of samples
First make sure that the dependencies are installed and then generate an ROI summary for each bam file.

### Creating a region of interest (ROI) summary for each bam file:
```Rscript ~/git_clone/CNV_Radar/bam2roi.r -b <bam file>.bam -d <bed file>.bed -z >> bam2roi.log 2>&1```

### Create a control dataset using normal samples:
```Rscript ~/git_clone/CNV_Radar/CNV_Radar_create_control.r directory=/Projects/SMM/normals/ out=/Projects/SMM/normals/CNV_SSV6_control.RData```

### Run CNV Radar:
```Rscript ~/git_clone/CNV_Radar/CNV_Radar.r -C /Projects/SMM/normals/CNV_SSV6_control.RData -S /Projects/SMM/Tumors/Sample1 -R _roiSummary.txt  -mmrfReport F -O /Projects/SMM/Tumors --includeScore T --highSensitivity T```

### Optional creation of ROI dendrograms:
```Rscript ~/git_clone/CNV_Radar/CreateROI_dendrograms.r bdir=/Projects/SMM/normals/ ds=SMM_normal```

## Outputs
There are two primary outputs from the tool: genomic plots in jpeg format and a <--out>.tsv output

|     Values     | Description                                                                                           |
|:--------------:|-------------------------------------------------------------------------------------------------------|
|       Chr      | Chromosome where the event occurs                                                                     |
|      Start     | Left most base pair in the CNV event                                                                  |
|      Stop      | Right most base pair in the CNV event                                                                 |
|     log2FC     | The log2 transformed fold change for the CNV event                                                    |
|     QScore     | The score used to identify CNV events [Bounded smoothed allele frequency x (smoothed fold change x 20)] |
| Observed Depth | The observed mean depth across the CNV event                                                          |
| Expected Depth | The expected mean depth across the CNV event                                                          |
|     Zscore     | (expected fold change - observed fold change) / standard deviation of fold change                     |
|     HetVar     | Mean heterozygous allele frequency across the CNV event                                               |
|      IsCNV     | T/F whether the CNV event is a significant                                                            |
|     Length     | Total number of base pairs in the CNV event                                                           |
|   IsLOH_Only   | T/F if the CNV event is only a copy number neutral loss of heterozygousity                            |

## Feedback
Yes! Please give us your feedback, raise issues, and let us know how the tool is working for you. Pull requests are welcome.
