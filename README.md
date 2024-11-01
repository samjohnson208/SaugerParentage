# SaugerParentage
Created by Sam Johnson, 10/25/24

sjohn208@uwyo.edu

## Project Summary

Starting the bioinformatics workflow for the Wind River Sauger Project! Whoop!

## Workflow Outline

* Demultiplex

* Split .fastq

* Alignment

* Variant Calling

* Filtering

* Entropy

## Demultiplex

### Unzip raw .fastq's

```{bash}
sbatch slurm_zip.sh
```

### Count raw reads

```{bash}
salloc --account=ysctrout --time=3:00:00
grep -c "^@" 1SaugEvens.fastq
   # 1,111,224,597
grep -c "^@" 1SaugOdds.fastq
   # 1,156,307,046
```

### Parsing

```{bash}
sbatch slurm_parse.sh
```