#!/bin/bash
#************************************************************************#
#BWA-MEM_Alignment Script by Lindsey Fenderson                           #
#Current version 1.3.1 - October 18, 2021                                #
#Pipeline for mapping fastq reads to a reference genome with the BWA-MEM # 
#  algorithm and additional quality control steps based on the GATK Best # 
#  Practices Data Pre-Processing for Variant Discovery Pipeline          #
#  (including samtools fixmate, index & sort; GATK Indel Realigner,      #
#  Picard Mark Duplicates, and GATK Base Quality Score Recalibration)    #
#  to generate analysis-ready bam files.                                 #
#************************************************************************#

#Program assumes that input data have already been quality processed (e.g., with the DataPreProcessingScript script) to evaluate the sequencing output and trim adapters and low quality bases or reads and that the files have been bzip2 compressed.
#To use the current version of this program, set up a parameter file with the following (see BWA-MEM_Alignment-Paramfile.txt for example): 
#   1.) Set the data path and directory of the files you wish to map. 
#       1A.) If samples were not preprocessed with the DataPreProcessingScript.sh, change the suffix of the batch of files to be mapped as needed in lines 31, 38 & 39.
#   2.) Change the path and reference genome name if required to the species you wish to map the reads to in line 56. Note that the reference genome must be bwa indexed (if it doesn't already exist, this can be generated in the linuxbrew/colsa module via e.g.,: bwa index -p [prefix] input_fasta_reference.fna. Example: bwa index -p Zalbicollis-1.0.1 GCF_000385455.1_Zonotrichia_albicollis-1.0.1_genomic.fna. The prefix is what the final pointer in the path should equal (i.e., don't have to provide entire reference genome name, just the prefix that was used to build the index)
#   3.) Change the output path if required to the corresponding directory of the species the reads have been mapped to in line 25. 
#   4.) Change the reference genome suffix that will be appended to the output filenames to the genome version the reads were mapped to if required in line 26.

#TO RUN -> Before running the script you need to 1.) assign the location of your parameter file as a variable, 2.) then add that argument to your sbatch job submission command. For example, you would use something like the following 2 command lines:
# export ParamFile=/mnt/lustre/mel/shared/Scripts/BWA-MEM_Alignment-Paramfile.txt
# sbatch ./LEF-BWA-MEM_AlignmentSlurmScript.sh –export=ParamFile
# ***It is recommended to use an informative slurm script which records additional pertinant details about the analysis run and importantly, names the slurm output informatively so you're not just left with a million unidentifiable slurm-#####.out files at the end of your project. See an example template at https://github.com/KovachLab/Bioinfo-Scripts/blob/main/BioinformaticBookkeepingSlurmTemplate.sh and the one I actually use at https://github.com/KovachLab/Bioinfo-Scripts/blob/main/LEF-FinalDataPreProcessingSlurmScript.sh. Snakemake is another good option for ensuring your analyses are well documented and reproducible.

module purge
module load linuxbrew/colsa
source ${ParamFile} 

#Map reads to reference genome
if [ $ReferenceGenome == "GECOSparrowSpecies" ]; then
    cd /mnt/lustre/mel/shared/GECO/Data/MappedData
    #Create a list of files to be processed
    Directory=$(echo $InputDataPath | rev | cut -d'/' -f 2 | rev)
    ls $InputDataPath*$QualityTrimmedReadSuffix > FilesToProcess$Directory
    sed -i 's,'"$InputDataPath"',,g' FilesToProcess$Directory

    #Create a root name list
    cut -d"." -f1 FilesToProcess$Directory > RootName$Directory
    sort RootName$Directory | uniq > UniqRootName$Directory
     echo -e "Mapping samples to their respective GECO species genomes (saltmarsh, Nelson's, seaside, song, swamp and/or savannah sparrow)"
     grep caudacuta UniqRootName$Directory > SALSUniqRootName$Directory
     grep nelsoni UniqRootName$Directory > NESPUniqRootName$Directory
     grep maritima UniqRootName$Directory > SESPUniqRootName$Directory
     grep melodia UniqRootName$Directory > SOSPUniqRootName$Directory
     grep georgiana UniqRootName$Directory > SWSPUniqRootName$Directory
     grep sandwichensis UniqRootName$Directory > SAVSUniqRootName$Directory
     for value in SALS NESP SESP SOSP SWSP SAVS
          do
          echo -e "Checking if $value samples exist."
          if [ `wc -l $value"UniqRootName"$Directory | awk '{print $1}'` -ge "1" ]; then
               echo -e "Mapping $value samples to $value genome."
               RefGenome=$(grep -w "$value"ReferenceGenome"" /mnt/lustre/mel/shared/Scripts/ReferenceGenomesList | cut -d"=" -f2); ReferenceGenome=$(sed -e 's/^"//' -e 's/"$//' <<<"$RefGenome")
                echo -e "Target reference genome is $ReferenceGenome"
                RefGenomeSuffix=$(grep -w "$value"ReferenceGenomeSuffix"" /mnt/lustre/mel/shared/Scripts/ReferenceGenomesList | cut -d"=" -f2); ReferenceGenomeSuffix=$(sed -e 's/^"//' -e 's/"$//' <<<"$RefGenomeSuffix")
                echo -e "Reference genome suffix is $ReferenceGenomeSuffix"
                OutDataPath=$(grep -w "$value"OutputDataPath"" /mnt/lustre/mel/shared/Scripts/ReferenceGenomesList | cut -d"=" -f2); OutputDataPath=$(sed -e 's/^"//' -e 's/"$//' <<<"$OutDataPath")
                echo -e "Output data path is $OutputDataPath"
                echo -e "Checking if samples already exist."
                ls ${OutputDataPath}*sorted.bam > ExistingMappedReads$Directory
                sed -i 's,'"$OutputDataPath"',,g' ExistingMappedReads$Directory
                sed -i 's,'"_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam"',,g' ExistingMappedReads$Directory
                cat $value"UniqRootName"$Directory ExistingMappedReads$Directory > FutureMappedSamples$Directory
                sort FutureMappedSamples$Directory | uniq -d > DuplicateLibraries$Directory
                comm -23 $value"UniqRootName"$Directory DuplicateLibraries$Directory > UniqLibraries$Directory

                if [ `wc -l "DuplicateLibraries" | awk '{print $1}'` -ge "1" ]; then
                    echo -e "Duplicate samples already exist. Assuming the duplicates are from separate libraries or sequencing runs and will merge mapped bams from the same individual"
                    mkdir -p $OutputDataPath/temp
        
                    while read DuplicateLibraries; do
                        File1="$DuplicateLibraries$QualityTrimmedReadSuffix1"
                        File2="$DuplicateLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        #Extract header information from fastq file
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        #Run bwa-mem
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${DuplicateLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                         #"Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags" #This command runs samtools fixmate, which fills in mate coordinates, ISIZE and mate related flags from a name-sorted alignment, and also sets the output format to bam. I opted to just use the defaults and did not include the flags to Remove secondary and unmapped reads (-r), Disable FR proper pair check (-p), Add template cigar ct tag (-c), Add ms (mate score) tags (-m) (These are used by markdup to select the best reads to keep), or Do not add a @PG line to the header of the output file (--no-PG). I suppose I would use the mate score tags if I used the samtools markdup tool instead of Picardtools MarkDuplicates...
                        samtools fixmate -m -O bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Merge duplicate bam files, re-sort the merged bam file and index the new merged bam file.
                        samtools merge ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        samtools sort -O bam -o ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.sorted.bam ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam
                        samtools index ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        ##rm ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam 
                        ##rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        rm ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam
                    done< DuplicateLibraries$Directory
                    
                    while read UniqLibraries; do
                        File1="$UniqLibraries$QualityTrimmedReadSuffix1"
                        File2="$UniqLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        #Extract header information from fastq file
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        #Run bwa-mem
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${UniqLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                         #"Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags" #This command runs samtools fixmate, which fills in mate coordinates, ISIZE and mate related flags from a name-sorted alignment, and also sets the output format to bam. I opted to just use the defaults and did not include the flags to Remove secondary and unmapped reads (-r), Disable FR proper pair check (-p), Add template cigar ct tag (-c), Add ms (mate score) tags (-m) (These are used by markdup to select the best reads to keep), or Do not add a @PG line to the header of the output file (--no-PG). I suppose I would use the mate score tags if I used the samtools markdup tool instead of Picardtools MarkDuplicates...
                        samtools fixmate -m -O bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                    done< UniqLibraries$Directory
                else
                    echo -e "No duplicate samples found. Continuing mapping pipeline."
                    while read UniqLibraries; do
                        File1="$UniqLibraries$QualityTrimmedReadSuffix1"
                        File2="$UniqLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        #Extract header information from fastq file
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        #Run bwa-mem
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${UniqLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                        #"Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags" #This command runs samtools fixmate, which fills in mate coordinates, ISIZE and mate related flags from a name-sorted alignment, and also sets the output format to bam. I opted to just use the defaults and did not include the flags to Remove secondary and unmapped reads (-r), Disable FR proper pair check (-p), Add template cigar ct tag (-c), or Do not add a @PG line to the header of the output file (--no-PG); in large part because I don't see these flags used, nor was I able to discern why employing these flags would be useful. Add ms (mate score) tags (-m) (These are used by markdup to select the best reads to keep),I suppose I would use the mate score tags if I used the samtools markdup tool instead of Picardtools MarkDuplicates...
                        samtools fixmate -m -O bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                    done< UniqLibraries$Directory
                fi
#Map input files to specified reference. Parameters include using 24 threads (or however many CPUs were requested in the slurm job script), the read group header line for the output SAM file (includes the sequencing instrument ID:run number on the instrument:flow cell ID from the header of the fastq files as the unique read group identifier; the root sample name of the sample that was sequenced; an identifier signifying the library batch the sequence was derived from; and the sequencing platform that generated this data), including alignments for single-end/unpaired paired-end reads so as not to discard any data, (but they will be flagged as secondary alignments), (want to use flag -C here to appending the barcode to the SAM file but is formatted incorrectly so as to lead to incorrect SAM output; need to figure out how to fix and e.g., make header in sam format with cs:Z: prefix; have removed flag for now), and mark shorter split hits as secondary for Picard tools compatibility. Then the prefix of the indexed genome is specified, followed by the R1 and R2 filenames to be mapped (specifying that they should be decompressed from bzip2 format to the standard output so the reads can be mapped by bwa) and the desired output filename. 

#Other parameters left as default (mostly because I don't have any understanding of if it is worthwhile or better to tune any of these parameters, and Heng Li is a pretty smart dude!) include: the minimum seed length of 19, band with of 100,Off-diagonal dropoff of 100, re-seeding value of 1.5, discard an alignment if it occurs more than 500 times in the genome, Drop chains shorter than 0.5 of the longest overlapping chain, perform at most 50 rounds of mate-SW, drop a chain if the number of bases in seeds is smaller than 0.  This option is primarily used for longer contigs/reads. When positive, it  also  affects seed filtering; did not use -P flag to rescue missing hits/ignore hits that fit a proper pair, default matching score of 1, default mismatch penalty of 4, default Gap open penalty of 6, default gap extension penalty of 1, clipping penalty of 5, unpaired read pair penalty of 17, did not use -p flag as input files are not interleaved, I didn't use the -q flag, which doesn't reduce the mapping quality of split alignment of lower alignment score, I didn't use the -5 flag which For split alignment, mark the segment with the smallest coordinate as the primary. It automatically applies option -q as well. This option  may  help  some Hi-C pipelines. By default, BWA-MEM marks highest scoring segment as primary; I didn't use the -T flag for not including alignment in output if it has a mapping quality score below 30 (this flag has no effect on paired end reads anyway apparently); I didn't use the -j flag which treats ALT contigs as part of the primary assembly (i.e. ignore the db.prefix.alt file); I didn't use the -h flag INT[,INT2] where if a query has not more than [5200] hits with score higher than 80% of the best hit, output them all in the XA tag. If INT2 is specified, BWA-MEM outputs up to INT2 hits if the list contains a hit to an ALT contig. [5,200]; I left hard clipping in place for the supplementary mappings, verbosity is what it is (all normal messages)!, and I didn't specify -I FLOAT[,FLOAT[,INT[,INT]], which specifies the mean, standard deviation (10% of the mean if absent), max (4 sigma from the mean if absent) and min (4 sigma if absent) of the insert size distribution. Only applicable to the FR orientation. By default, BWA-MEM infers these numbers and the pair orientations given enough reads. [inferred]
      # Filter aligned reads with a mapping quality (Phred) below this value
##      MinQuality: 0
      # Filter reads that did not map to the reference sequence
##      FilterUnmappedReads: yes

          else
               echo -e "No $value samples found."
          fi
     done
 #Clean up temp files
     ##rm RootName
     ##rm *UniqRootName
     ##rm DuplicateLibraries UniqLibraries
     ##rm ExistingMappedReads FutureMappedSamples
else
     echo -e "Mapping samples to $ReferenceGenomeSuffix genome."
#Map reads to reference genome
cd $OutputDataPath
#Create a list of files to be processed
    Directory=$(echo $InputDataPath | rev | cut -d'/' -f 2 | rev)
    ls $InputDataPath*$QualityTrimmedReadSuffix > FilesToProcess$Directory
    sed -i 's,'"$InputDataPath"',,g' FilesToProcess$Directory

    #Create a root name list
    cut -d"." -f1 FilesToProcess$Directory > RootName$Directory
    sort RootName$Directory | uniq > UniqRootName$Directory
                echo -e "Target reference genome is $ReferenceGenome"
                echo -e "Output data path is $OutputDataPath"
                echo -e "Checking if samples already exist."
                ls ${OutputDataPath}*sorted.bam > ExistingMappedReads$Directory
                sed -i 's,'"$OutputDataPath"',,g' ExistingMappedReads$Directory
                sed -i 's,'"_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam"',,g' ExistingMappedReads$Directory
                cat UniqRootName$Directory ExistingMappedReads$Directory > FutureMappedSamples$Directory
                sort FutureMappedSamples$Directory | uniq -d > DuplicateLibraries$Directory
                comm -23 UniqRootName$Directory DuplicateLibraries$Directory > UniqLibraries$Directory
                
                if [ `wc -l "DuplicateLibraries"$Directory | awk '{print $1}'` -ge "1" ]; then
                    echo -e "Duplicate samples already exist. Assuming the duplicates are from separate libraries or sequencing runs and will merge mapped bams from the same individual"
                    mkdir -p $OutputDataPath/temp
        
                    while read DuplicateLibraries; do
                        File1="$DuplicateLibraries$QualityTrimmedReadSuffix1"
                        File2="$DuplicateLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${DuplicateLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                        samtools fixmate -m -O bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Merge duplicate bam files, re-sort the merged bam file and index the new merged bam file.
                        samtools merge ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        samtools sort -O bam -o ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.sorted.bam ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam
                        samtools index ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        ##rm ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam 
                        ##rm ${OutputDataPath}"temp/"${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        rm ${OutputDataPath}${DuplicateLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.merged.bam
                    done< DuplicateLibraries$Directory

                    while read UniqLibraries; do
                        File1="$UniqLibraries$QualityTrimmedReadSuffix1"
                        File2="$UniqLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${UniqLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                        samtools fixmate -m -O bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                    done< UniqLibraries$Directory
                else
                    echo -e "No duplicate samples found. Continuing mapping pipeline."
                    while read UniqLibraries; do
                        File1="$UniqLibraries$QualityTrimmedReadSuffix1"
                        File2="$UniqLibraries$QualityTrimmedReadSuffix2"
                        echo -e "The following files are being mapped: File1 = $File1 & File2 = $File2"
                        Read1="$InputDataPath$File1"
                        echo -e "Files being mapped with their full paths:"
                        echo -e "Read1 = $Read1"
                        Read2="$InputDataPath$File2"
                        echo -e "Read2 = $Read2"
                        #Extract header information from fastq file
                        header=$(less $Read1 | head -n 1)
                        id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/M_//g')
                        bwa mem -t $SLURM_CPUS_PER_TASK -R "@RG\tID:${id}\tSM:${UniqLibraries}\tLB:${Library}\tPL:${Platform}" -a -q -Y -K 80000000 $ReferenceGenome <(bzip2 -dc $Read1) <(bzip2 -dc $Read2) > "${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam"
                        #"Because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags" #This command runs samtools fixmate, which fills in mate coordinates, ISIZE and mate related flags from a name-sorted alignment, and also sets the output format to bam. I opted to just use the defaults and did not include the flags to Remove secondary and unmapped reads (-r), Disable FR proper pair check (-p), Add template cigar ct tag (-c), Add ms (mate score) tags (-m) (These are used by markdup to select the best reads to keep), or Do not add a @PG line to the header of the output file (--no-PG). I suppose I would use the mate score tags if I used the samtools markdup tool instead of Picardtools MarkDuplicates...
                        samtools fixmate -m -O bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        #Most subsequent tools require coordinate-sorted (instead of read name-sorted) bam files and bam files need to be indexed to facilitate their use with various tools, hence the next 2 commands.
                        samtools sort -O bam -o ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                        samtools index ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam
                        #Cleanup temp files
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.sam
                        rm ${OutputDataPath}${UniqLibraries}_trimmed-bwamem-$ReferenceGenomeSuffix.bam
                    done< UniqLibraries$Directory
                fi
     #Clean up temp files
     ##rm RootName
     ##rm UniqRootName
     ##rm DuplicateLibraries UniqLibraries
     ##rm ExistingMappedReads FutureMappedSamples
fi

#Run GATK Indel Realigner
#Note that this tool requires an indexed reference genome (if it doesn't already exist, this can be generated in the linuxbrew/colsa module via e.g.,: samtools faidx GCF_000002315.6_GRCg6a_genomic.fna) and a reference sequence dictionary (if it doesn't already exist, this can be generated in the linuxbrew/colsa module via e.g.,: gatk CreateSequenceDictionary -R GCF_000002315.6_GRCg6a_genomic.fna)
#Indel realigner is not supported in GATK4, as from what I can gather, the tool has been incorporated into the GATK vcf callers. I decided to include this step regardless to better permit the use of alternate variant discovery tools that don't necessarily incorporate realignment. Thus, to utilize the tool I need to switch to a conda environment so can use GATK3.8.
module purge
module load anaconda/colsa

#Use the RealignerTargetCreator tool to generate a set of sites likely to contain indels
##while read UniqRootName; do
##java -jar /mnt/lustre/mel/shared/Software/GenomeAnalysisTK.jar \
##    -T RealignerTargetCreator \
##    -R /mnt/lustre/mel/shared/GECO/ReferenceGenomes/DomesticChicken/GCF_000002315.6_GRCg6a_genomic.fna \
##    -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam \
##    -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.intervals
#Perform local realignment around indels
##java -Xmx8G -Djava.io.tmpdir=/tmp -jar /mnt/lustre/mel/shared/Software/GenomeAnalysisTK.jar \
##    -T IndelRealigner --filter_bases_not_stored \
##    -R /mnt/lustre/mel/shared/GECO/ReferenceGenomes/DomesticChicken/GCF_000002315.6_GRCg6a_genomic.fna \
##    -targetIntervals ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.intervals \
##    -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.bam \
##    -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.bam
##done< UniqRootName

#Mark PCR duplicates
module purge
module load linuxbrew/colsa
#***Note that using Picardtools probably under represents the number of PCR duplicates since I am using collapsed reads. In the future should probably incorporate the paleomix rmdup_collapsed tool (This tool "Filters PCR duplicates for merged/collapsed paired-ended reads, such as those generated by AdapterRemoval with the --collapse option enabled. Unlike SAMtools rmdup or Picard MarkDuplicates, this tool identifies duplicates based on both the 5' and the 3' alignment coordinates of individual reads.") - paleomix rmdup_collapsed --remove-duplicates < sorted.bam > < out.bam >
##while read UniqRootName; do
##gatk MarkDuplicates -I ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.bam -O ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.marked.indelrealigned.bam -M ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.indelrealigned.dup_metrics.txt
#Resort the files
#I'm unclear as to why this needs to be done again, but various pipelines I've seen online do it, so I'm following the herd off the cliff...can't hurt at any rate!
##samtools sort -O bam -o ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.resorted.marked.indelrealigned.bam ${UniqRootName}_trimmed-bwamem-$ReferenceGenomeSuffix.sorted.marked.indelrealigned.bam
##done< UniqRootName

#At this stage it is ideally recommended to rescale the base quality scores. This requires a known variants file. I used a chicken vcf from http://bigd.big.ac.cn/chickensd/, which was developed from 863 chicken genomes (reported in Wang et al. (2020) 863 genomes reveal the origin and domestication of chicken. Cell Research), as it makes more sense to me to use a chicken reference in this case since all 6 species were mapped to the chicken reference genome, than to use the only sparrow vcfs known to me (from Jen's 2019 EcoEvo paper for STSP only and 1 SESP, and from Jen's 2019 EvolLett paper for SAVS, NESP, SOSP and SWSP). As I understand it, this step can also be omitted.


#Cleanup temp files
##/mnt/lustre/mel/leq29/BioinformaticScripts/rm *.intervals
##/mnt/lustre/mel/leq29/BioinformaticScripts/rm *.sorted.bam
##/mnt/lustre/mel/leq29/BioinformaticScripts/rm *.sorted.marked.indelrealigned.bam
##/mnt/lustre/mel/leq29/BioinformaticScripts/rm *.sorted.indelrealigned.bam
##/mnt/lustre/mel/leq29/BioinformaticScripts/rm *.sorted.bam.bai
#Also need to fix the read group/library stuff so it extracts the info from a Yaml 

