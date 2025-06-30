#!/bin/bash

# this script is by MPR. working to edit my own. added to git 063025.
# calculates the number of reads, number of loci, mean sequencing depth and number of loci with min 10 reads and their mean depth  


#SBATCH --job-name=pcrdup
#SBATCH -A wagnerlab
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10   
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mrodri23@uwyo.edu


module spider | head
module load arcc/1.0 gcc/14.2.0 samtools/1.20

path=$1
lib=$2

cd "${path}"

# create the report.txt file with a line for each individual (each bam file) with read counts
echo -e "sample\tsampleLib\tmappedReads" > report.txt

for i in $(ls *\.sorted.bam);
do
echo "counting reads in "$i
totalReads=`samtools view $i | wc -l`
echo -e ${i%.sorted.bam}"\t"${lib}"\t"$totalReads >> report.txt;
done

# creates a file with a line for each locus with mapped reads in each contig for each individual
echo "sample marker scaffold locus orientation depth" > seq_depth.txt

for i in *.sorted.bam; do
    echo "Working on $i"
    sample=$(basename "$i" .sorted.bam)

    samtools view "$i" | awk -v sample="$sample" '
    BEGIN {
        markerold = ""
        marker = ""
        scaffold = ""
        locus = ""
        orientation = ""
        depth = ""
    }
    {
        if ($5 > 30) {
            if (scaffold == "") {
                # First read on first scaffold
                if ($2 == 16) {
                    b = 0
                    cigar = $6
                    gsub(/M/, "+", cigar)
                    gsub(/D/, "+", cigar)
                    gsub(/[0-9]+I/, "", cigar)
                    sub(/$/, "0", cigar)
                    split(cigar, indexx, "+")
                    for (i in indexx) b += indexx[i]
                    locus = $4 + b - 4
                } else {
                    locus = $4
                }
                marker = $3 "_" locus "_" $2
                scaffold = $3
                orientation = $2
                depth = 1
            }
            else if (scaffold == $3) {
                markerold = marker
                locusold = locus

                if ($2 == 16) {
                    b = 0
                    cigar = $6
                    gsub(/M/, "+", cigar)
                    gsub(/D/, "+", cigar)
                    gsub(/[0-9]+I/, "", cigar)
                    sub(/$/, "0", cigar)
                    split(cigar, indexx, "+")
                    for (i in indexx) b += indexx[i]
                    locus = $4 + b - 4
                } else {
                    locus = $4
                }
                marker = $3 "_" locus "_" $2

                if (marker == markerold) {
                    depth += 1
                } else {
                    print sample, markerold, scaffold, locusold, orientation, depth
                    scaffold = $3
                    if ($2 == 16) {
                        b = 0
                        cigar = $6
                        gsub(/M/, "+", cigar)
                        gsub(/D/, "+", cigar)
                        gsub(/[0-9]+I/, "", cigar)
                        sub(/$/, "0", cigar)
                        split(cigar, indexx, "+")
                        for (i in indexx) b += indexx[i]
                        locus = $4 + b - 4
                    } else {
                        locus = $4
                    }
                    marker = $3 "_" locus "_" $2
                    orientation = $2
                    depth = 1
                }
            } else {
                print sample, marker, scaffold, locus, orientation, depth
                scaffold = $3
                if ($2 == 16) {
                    b = 0
                    cigar = $6
                    gsub(/M/, "+", cigar)
                    gsub(/D/, "+", cigar)
                    gsub(/[0-9]+I/, "", cigar)
                    sub(/$/, "0", cigar)
                    split(cigar, indexx, "+")
                    for (i in indexx) b += indexx[i]
                    locus = $4 + b - 4
                } else {
                    locus = $4
                }
                marker = $3 "_" locus "_" $2
                orientation = $2
                depth = 1
            }
        }
    }
    END {
        print sample, marker, scaffold, locus, orientation, depth
    }' >> seq_depth.txt;
done

awk '{if($6>9) print $0}' seq_depth.txt > seq_depth_min10.txt

# get per locus information
echo -e "locus\ttotalInds\ttotalReads" > contigs.txt
cat seq_depth.txt | \
awk '{
 counter[$2]++
 totReads[$2]+=$6
}
END{
 for (i in counter) print i"\t"counter[i]"\t"totReads[i]
}' >> contigs.txt

# get per locus information (count only individuals with at least 10 reads)
echo -e "locus\ttotalInds\ttotalReads" > lociMin10reads.txt
cat seq_depth_min10.txt | \
awk '{
 counter[$2]++
 totReads[$2]+=$6
}
END{
 for (i in counter) print i"\t"counter[i]"\t"totReads[i]
}' >> lociMin10reads.txt

# create a file with lociN and meanDepth per ind
echo -e "sample lociN meanDepth" > indInfo.txt
awk '{
counter[$1]++
depth[$1]+=$6
}
END{
for(i in counter) print i,counter[i],depth[i]/counter[i]
}' seq_depth.txt >> indInfo.txt


# create a file with lociN and meanDepth per ind counting loci with min 10 reads
echo -e "sample lociN meanDepth" > indInfoMin10reads.txt
awk '{
counter[$1]++
depth[$1]+=$6
}
END{
for(i in counter) print i,counter[i],depth[i]/counter[i]
}' seq_depth_min10.txt >> indInfoMin10reads.txt

# sort the stats files skipping the header line (careful: "sample" cannot be included in any sample name)
sort indInfoMin10reads.txt | awk '!/sample/' > sorted_indInfoMin10reads.txt
sort indInfo.txt | awk '!/sample/' > indInfoSorted
sort report.txt | awk '!/sample/' > sortedReport


# create a file with Nreads, lociN, meanDepth... for each individual
join -a 1 -a 2 sortedReport indInfoSorted -1 1 -2 1 > indInfoAdv.txt
sort indInfoAdv.txt > indInfoAdvSorted
echo "sample sampleLib mappedReads lociN meanDepth lociNmin10reads meanDepthMin10reads" > mappingReport.txt
join -a 1 -a 2 indInfoAdvSorted sorted_indInfoMin10reads.txt -1 1 -2 1 >> mappingReport.txt


# Remove files not needed anymore:
rm indInfoAdvSorted indInfoAdv.txt sortedReport indInfoSorted sorted_indInfoMin10reads.txt contigs.txt lociMin10reads.txt indInfoMin10reads.txt indInfo.txt report.txt 
