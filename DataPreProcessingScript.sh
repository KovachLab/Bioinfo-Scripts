#!/bin/bash
#*********************************************************************#
#Data Preprocessing Script by Lindsey Fenderson                       #
#Current version 1.3.0 - May 17, 2021                                 #
#Script for creating fastqc reports on newly sequenced data,          #
#  (code for double checking barcodes if needed),                     #
#  & trimming adapters with AdapterRemoval2 on UNH's Premise cluster. #
#*********************************************************************#


#Process assumes that sequencing data for the same individual generated on multiple lanes in the same sequencing run has already been concatenated and that the raw sequencing files have been named informatively according to the Kovach_Lab_File_Naming_SOPs (i.e., expected input files have 4 fields separated by underscores before the R1/R2 suffix, e.g.: IndividualID_Species_SamplingLocation_CaptureDate_R#.f*q - if this is not the case, the script will likely not work as expected and you'll need to edit throughout e.g., if the field separator isn't an underscore, or if there are a different number of fields before the R1 or R2 suffix etc...

##Expected output of this program: 
# 1.) This program will generate fastqc analysis reports on the individual R1 and R2 raw sequencing data and output them to the folder 'fastqcRawData'.
# 2.) The program then quality trims adapters and low quality termini from reads and collapses overlapping read pairs using AdapterRemoval2. Output files will have the same root name as the input fastq files, with the suffix .pair#.truncated.bz2 (where # is the number 1 or 2 for reads 1 or 2 respectively). This script also generates charts and reports summarizing the read number and length distribution statistics before and after trimming, in the directory AdapterRemovalSettings.
# 3.) The program will then generate fastqc analysis reports on the individual R1 and R2 quality trimmed sequencing data generated in step 2 and output them to the folder 'fastqcTrimmedData'. 
# 4.) The program also runs multiQC and generates summary charts and reports for the bcl2fastq data (if available), the raw and trimmed fastQC reports, and the Adapter Removal files.
# 5.) Finally, the program also performs some data and file cleanup and compression to best preserve disc space; archiving all of the raw reads along with any sequencer metadata and the adapter removal settings files in $Directory_RawSequencingReads.tar.bz2, and archiving all of the multiQC and fastQC output for the raw and trimmed reads in $Directory_QCReports.tar.bz2. These files, along with all of the pdf Adapter Removal Settings files and QC reports are organized into the directory "RawReads-and-SequencingPreProcessingData". Copies of the pdf output are automatically backed up to the cloud (e.g., Box) and LaTeX code written to easily import and index the pdfs into an electronic lab notebook. The script also records the fastqc, multiqc and AdapterRemoval program version numbers that were used near the top of the slurm output file.

##To use this program:
#   1.) Program assumes the raw reads are in fastq format (files can be compressed (any format, e.g., gzip or bzip2) or not) with a suffix involving fq or fastq. If not, change the suffix as needed in lines 41 and 445.
#   2.)  FASTQC parameters should be provided in associated files at the paths listed in the paramfile. Parameters used will print to the slurm stdout file. Adapter Removal parameters should also be given in the paramfile; other parameters automatically used are described below in line 133; Adapter Removal paramaters used are recorded in the output .settings files and AdapterRemovalParameters pdfs.
#   3.) Edit "Reports" and "Stats" and/or add the directory names or add any other files as needed for the raw sequencing metadata to be archived with the raw reads as relevant for your data (e.g, "report", "statsPlate1", "MD5.txt", "Summary.txt", "checkSize.xls", "Rawdata_Readme.pdf" etc.) in lines 421-422(+). Be sure to include the absolute or relative path if for some reason the files or directories are not in the working directory.
#   4.) In your working directory, rename your undetermined/unknown fastq files to have the same naming structure as the other files to make the script happy (i.e., the name must start with "Undetermined" (or else edit the prefix in line 41), then add 3 additional dummy fields for the "species", "sampling location" and "capture date" fields, separated by an underscore, e.g.: Undetermined_UNK_NA_S0_R1.fastq.gz)


#   ***Note 1: As per our Kovach_Lab_Data_Management_User_Guide-START_HERE SOP, it is required that the raw sequencing reads are backed up in 2 separate locations (both in the cloud and on a physical hard drive located on campus) and that the integrity of the backups is confirmed PRIOR to running this script, since the combined raw read file archive may be too large to effectively upload or download to other backup locations (or it could get corrupted in the process) and the individual raw read files will be deleted once they are added to the archive (so you should still be able to access the raw reads if you decompress and unarchive the $Directory_RawSequencingReads.tar.bz2 file, however there is always a risk of data loss if the file gets corrupted). As noted above, the raw reads also need to be concatenated and named according to our SOPs.
source /mnt/lustre/mel/shared/Scripts/DataPreProcessingScriptv1.3.0-Paramfile.txt

module purge
module load linuxbrew/colsa

#Record program versions:
echo -e "Program versions used in this script: \n"
fastqc -v
AdapterRemoval --version
module purge
module load anaconda/colsa
conda activate multiqc-1.10.1
multiqc --version
module purge
module load linuxbrew/colsa

echo -e "\n"
##Generate FastQC reports for your sequencing data##
mkdir ./fastqcRawData
#Create a list of files to be processed (DOES process Undetermined raw fastq files for fastqc)
ls *.f*q* > AllFilesToProcess
awk '!/Undetermined*/' AllFilesToProcess > FilesToProcess 
#Calculate read length and max Ns
AnyFile=$(head -n 1 AllFilesToProcess)
ExampleRead=$(less $AnyFile |head -n 2 | tail -n 1 )
ReadLength=$(expr length $ExampleRead)
CycleLength=$(echo "($ReadLength - 1)" | bc)
Maxn=$(echo "($CycleLength * $PctNs)" | bc)
MaxNs=$( printf "%.0f" $Maxn )
#minalignmentlength=11
#Run fastqc
echo -e "\n *** RUNNING FASTQC ON RAW READS *** \n\n"
echo -e "FastQC parameter limits for raw read data:"
while read RawFastQCFile; do
    echo -e $RawFastQCFile
done< $RawFastQCLimitsFile  

while read AllFilesToProcess; do
    fastqc --limits $RawFastQCLimitsFile --adapters $AdapterFile --contaminants $ContaminantFile --kmers $kmersize -o fastqcRawData -t 24 -f fastq $AllFilesToProcess
done< AllFilesToProcess

#Create a root name list
cut -d"_" -f1-4 AllFilesToProcess > RootName
sort RootName | uniq > UniqRootName

line1=$(head -n 1 AllFilesToProcess)
Suffixa1=$(echo $line1 | cut -d"_" -f5)
Suffix1="_$Suffixa1"
line2=$(head -n 2 AllFilesToProcess | tail -n 1)
Suffixa2=$(echo $line2 | cut -d"_" -f5)
Suffix2="_$Suffixa2"
Suffix1A='_R1_fastqc.html'
Suffix2A='_R2_fastqc.html'
Suffix1B='_R1_fastqc.pdf'
Suffix2B='_R2_fastqc.pdf'

#Parse directory name & start LaTeX code for ELN
Directory=$(pwd | rev | cut -d'/' -f 1 | rev)
SeqDate=$(echo $Directory | cut -d'_' -f 1 )
SeqPlatform=$(echo $Directory | cut -d'_' -f 2 )
SeqPerson=$(echo $Directory | cut -d'_' -f 3 )
SeqDescription=$(echo $Directory | cut -d'_' -f 4 )
touch $Directory-ELNLaTeXCodeFastQCRawReads.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
echo '\begin{landscape} %Rotates the below pages to landscape orientation. Note also the 'angle=-270' command added to the \includepdf options; this inserts the pdf into landscape orientation on the page as well.' >> $Directory-ELNLaTeXCodeFastQCRawReads.txt

while read UniqRootName; do
    #Print all of the fastqc reports to pdf
    File1A="$UniqRootName$Suffix1A"
    File2A="$UniqRootName$Suffix2A"
    File1B="$UniqRootName$Suffix1B"
    File2B="$UniqRootName$Suffix2B"
    pandoc fastqcRawData/$File1A -V geometry:landscape -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o $File1B
    pandoc fastqcRawData/$File2A -V geometry:landscape -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o $File2B
    #Automatically generate LaTeX code to insert pdf fastqc reports into electronic lab notebook
    echo '\phantomsection' >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    #Parse sample name
    SampleID=$(echo $UniqRootName | cut -d'_' -f 1 )
    Species=$(echo $UniqRootName | cut -d'_' -f 2 )
    SampleLoc=$(echo $UniqRootName | cut -d'_' -f 3 )
    CaptureDate=$(echo $UniqRootName | cut -d'_' -f 4 )
    echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!FastQC Reports\!Raw data|(}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/QCReports/FastQCReports/fastqcRawData/$UniqRootName$Suffix1B}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/QCReports/FastQCReports/fastqcRawData/$UniqRootName$Suffix2B}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!FastQC Reports\!Raw data|)}" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
    echo "" >> $Directory-ELNLaTeXCodeFastQCRawReads.txt
done< UniqRootName
echo '\end{landscape} %Stops rotating pages to landscape orientation.' >> $Directory-ELNLaTeXCodeFastQCRawReads.txt

#Recreate root name list without the Undetermined files
cut -d"_" -f1-4 FilesToProcess > RootName
sort RootName | uniq > UniqRootName

echo -e "\n\n *** RUNNING ADAPTER-REMOVAL ON RAW READS *** \n\n"

#Run AdapterRemoval2 (does not process Undetermined fastq files)
#Echo to stdout record AdapterRemoval parameters not recorded in settings file used in analysis:
echo ""
echo -e "ADAPTER REMOVAL Maximum Number of Ns/Read Permitted Parameter = $MaxNs"
echo ""

while read UniqRootName; do
     File1="$UniqRootName$Suffix1"
     File2="$UniqRootName$Suffix2"
     #Optional - Double check the adapters before trimming
     #If needed, you can confirm the adaptors to search for and remove with the following command (I haven't integrated it into the script yet, but you can just run manually if you are unsure of the adapter sequence used; then just replace the consensus output with the adapters used below in the next step):
     ####AdapterRemoval --identify-adapters --file1 $File1 --file2 $File2
     #FYI the adapters used below are standard for Illumina Truseq dual-indexed libraries (adapter 1 written as 5'-3'; adapter 2 always given as the reverse complement sequence), containing the Illumina Universal adapter sequencing primer binding site sequence (read 1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC; read 2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT), followed by the sample's unique 8-mer index, and the i7/i5 primer sequence in every Illumina adapter that binds to the flow cell (adapter 1: CGTATGCCGTCTTCTGCTTG; adapter 2: TCGGTGGTCGCCGTATCATT). Note that this also removes the read 1 3' A overhangs added during library prep (hence why there is an 'A' at the start of the read 1 adapter sequence before the primer binding site sequence). One of the other most common Illumina Truseq adapter combinations used is just a single index adapter, in which case just delete 'NNNNNNNNGTGT' from the middle of adapter 2; also note that the length of the indexing barcode has changed over time so your adapter 1 may have fewer N's (e.g., 6 (rather than the 8 used here) is also common).
     
#Trim Adapters. Parameters extracted from the parameter file include the minimum base call quality threshold, the size of the window used for assessing average base call quality, the minimum required read length below which the read is discarded, and the minimum number of required basepair overlap in order to merge read 1 and read 2 pairs. Other parameters used include setting the program to use 24 threads, the fastq read number in the input files is separated from the header by a space, trim consecutive Ns from the 5' and 3' termini, and trim consecutive stretches of low quality bases at the 5' and 3' termini (below the set $MinQual threshold), discard reads that have X% (set in param file, used to calculate $MaxNs which is dependent on read length) or more Ns, and require reads to be at least $MinLength bp after trimming. The adapter sequences to be trimmed from the 3' ends of the reads are given, allowing ambiguity for the individual sample index so the code will run across samples. Read pairs that overlap by a minimum of $minalignmentlength bases are collapsed into a single contig. Output files are combined, meaning all of the trimmed read1 and read2 reads as well as collapsed read1/read2 reads (added to the read1 file) and singletons where one of the mate pairs was discarded (with reads discarded due to quality filters or read merging are replaced with a single 'N' and Phred score of 0) are included in the single output read 1 and read 2 files. These output files are  bzip2 compressed at the max compression level (9). Remaining parameters were left as default (i.e., --mm mismatch rate of 3, --shift of 2 bp in read2 to allow for missing bases, essentially no limit on max read length, and the default --minalignmentlength of 11bp overlap required before collapsing read1 and read2 pairs.)

     AdapterRemoval --file1 $File1 --file2 $File2 --basename $UniqRootName --threads 24 --mate-separator ' ' --trimwindows $TrimWin --trimns --trimqualities --minquality $MinQual --maxns $MaxNs --minlength $MinLength --minalignmentlength $minalignmentlength --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --collapse --combined-output --bzip2
done< UniqRootName

##Generate FastQC reports for your trimmed sequencing data##
mkdir fastqcTrimmedData

#Create a list of files to be processed (does not process Undetermined fastq files)
ls --ignore=Undetermined* | grep .truncated.bz2 > TrimmedFiles
echo -e "\n\n *** RUNNING FASTQC ON TRIMMED READS *** \n\n"

echo -e "FastQC parameter limits for quality trimmed read data:"
while read TrimmedFastQCFile; do
    echo -e $TrimmedFastQCFile
done< $QTrimmedFastQCLimitsFile  

while read TrimmedFiles; do
     fastqc --limits $QTrimmedFastQCLimitsFile --adapters $AdapterFile --contaminants $ContaminantFile --kmers $kmersize -o fastqcTrimmedData -t 24 -f fastq $TrimmedFiles
done< TrimmedFiles

#Start LaTeX code for ELN
touch $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
echo '\begin{landscape} %Rotates the below pages to landscape orientation. Note also the 'angle=-270' command added to the \includepdf options; this inserts the pdf into landscape orientation on the page as well.' >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
touch $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt

#Print the fastqc reports & adapter removal settings to pdf
while read UniqRootName; do
    Suffix1C='.pair1.truncated_fastqc.html'
    Suffix2C='.pair2.truncated_fastqc.html'
    Suffix1D='.pair1.truncated_fastqc.pdf'
    Suffix2D='.pair2.truncated_fastqc.pdf'
    Suffix1E='.settings'
    Suffix1F='.AdapterRemovalSettings.pdf'
    File1C="$UniqRootName$Suffix1C"
    File2C="$UniqRootName$Suffix2C"
    File1D="$UniqRootName$Suffix1D"
    File2D="$UniqRootName$Suffix2D"
    File1E="$UniqRootName$Suffix1E"
    File1F="$UniqRootName$Suffix1F"
    pandoc fastqcTrimmedData/$File1C -V geometry:landscape -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o $File1D
    pandoc fastqcTrimmedData/$File2C -V geometry:landscape -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o $File2D
    
    #Automatically generate LaTeX code to insert pdf fastqc reports into electronic lab notebook
    echo '\phantomsection' >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    #Parse sample name
    SampleID=$(echo $UniqRootName | cut -d'_' -f 1 )
    Species=$(echo $UniqRootName | cut -d'_' -f 2 )
    SampleLoc=$(echo $UniqRootName | cut -d'_' -f 3 )
    CaptureDate=$(echo $UniqRootName | cut -d'_' -f 4 )
    echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!FastQC Reports\!Trimmed data|(}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "\phantomsection" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!Adapter \& Quality Trimming Reports|(}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/QCReports/FastQCReports/fastqcTrimmedData/$UniqRootName$Suffix1D}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/QCReports/FastQCReports/fastqcTrimmedData/$UniqRootName$Suffix2D}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "\includepdf[pages=-,frame,scale=0.85,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/AdapterRemovalSettings/$UniqRootName-AdapterRemovalParameters.pdf}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo '\begin{landscape}' >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/AdapterRemovalSettings/$UniqRootName.AdapterRemovalDistributionCharts.pdf}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo '\end{landscape}' >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!FastQC Reports\!Trimmed data|)}" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "\index[sample]{$SampleID\_$Species\_$SampleLoc\_$CaptureDate\!Adapter \& Quality Trimming Reports|)}" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
    echo "" >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt
    echo "" >> $Directory-ELNLaTeXCodeAdapterRemovalSettings.txt
done< UniqRootName
echo '\end{landscape} %Stops rotating pages to landscape orientation.' >> $Directory-ELNLaTeXCodeFastQCTrimmedReads.txt

#Make Adapter Removal reports
LinesOut=$(echo "($minalignmentlength-2)" | bc)
PEReadLength=$(echo "($ReadLength * 2)-$LinesOut" | bc)
while read UniqRootName; do
    Suffix1E='.settings'
    File1E="$UniqRootName$Suffix1E"
    head -n 46 $File1E > $File1E-parameters
    tail -n $PEReadLength $File1E > $File1E-R.txt
done< UniqRootName
echo -e "\n\n *** RUNNING AdapterRemovalSettingsCharts.r SCRIPT IN R *** \n\n"
ls *R.txt > RInputFiles
while read RInputFiles; do
    FILENAME=${RInputFiles}
    Rscript /mnt/lustre/mel/shared/Scripts/AdapterRemovalSettingsCharts.r ${FILENAME} --save
done< RInputFiles

while read UniqRootName; do
    echo "<html><body><h1><b>$UniqRootName</b></h1></body></html>" >> $UniqRootName-AdapterRemovalParameters
    head -n 2 $UniqRootName.settings-parameters >> $UniqRootName-AdapterRemovalParameters
    echo '<br />' >> $UniqRootName-AdapterRemovalParameters
    head -n 6 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        A1="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>ADAPTER SEQUENCES</b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo ""  >> $UniqRootName-AdapterRemovalParameters
        echo '<html><body><p><b>Adapter1[1]: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$A1" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 7 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        A2="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Adapter2[1]: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$A2" >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<html><body><p><b>ADAPTER TRIMMING</b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo ""  >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 11 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        RNG="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>RNG seed: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$RNG" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 12 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        shift="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Alignment shift value: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$shift" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 13 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        mismatch="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Global mismatch threshold: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$mismatch" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 14 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Qualf="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Quality format (input): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Qualf" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 15 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Qualsi="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Quality score max (input): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Qualsi" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 16 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Qualfo="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Quality format (output): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Qualfo" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 17 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Qualso="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Quality score max (output): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Qualso" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 18 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Matesep="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Mate-number separator (input): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Matesep" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 19 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Trim5p="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Trimming 5p: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Trim5p" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 20 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Trim3p="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Trimming 3p: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Trim3p" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 21 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        TrimN="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Trimming Ns: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$TrimN" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 22 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        TrimPhred="$(echo $AdapterRemoval | cut -d'=' -f2)"
        echo '<html><body><p><b>Trimming Phred scores <= </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$TrimPhred" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 23 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        TrimWin="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Trimming using sliding windows: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$TrimWin" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 24 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        MinGL="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Minimum genomic length: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$MinGL" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 25 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        MxGL="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Maximum genomic length: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$MxGL" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 26 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Coll="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Collapse overlapping reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Coll" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 26 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Mino="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Minimum overlap (in case of collapse): </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Mino" >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<html><body><p><b>TRIMMING STATISTICS</b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo ""  >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 31 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        TotRP="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Total number of read pairs: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$TotRP" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 32 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        NuRP="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of unaligned read pairs: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$NuRP" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 33 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        NwaRP="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of well aligned read pairs: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$NwaRP" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 34 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Ndm1="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of discarded mate 1 reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Ndm1" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 35 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nsm1="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of singleton mate 1 reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nsm1" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 36 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Ndm2="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of discarded mate 2 reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Ndm2" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 37 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nsm2="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of singleton mate 2 reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nsm2" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 38 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nra1="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of reads with adapters[1]: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nra1" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 39 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nflcp="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of full-length collapsed pairs: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nflcp" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 40 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Ntcp="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of truncated collapsed pairs: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Ntcp" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 41 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nrr="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of retained reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nrr" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 42 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        Nrn="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Number of retained nucleotides: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$Nrn" >> $UniqRootName-AdapterRemovalParameters
    done
    head -n 43 $UniqRootName.settings-parameters | tail -n 1 | while read AdapterRemoval; do
        ALrr="$(echo $AdapterRemoval | cut -d':' -f2)"
        echo '<html><body><p><b>Average length of retained reads: </b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo -e "$ALrr" >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<br />' >> $UniqRootName-AdapterRemovalParameters
        echo '<html><body><p><b>LENGTH DISTRIBUTION CHARTS:</b></p></body></html>' >> $UniqRootName-AdapterRemovalParameters
        echo ""  >> $UniqRootName-AdapterRemovalParameters
        pandoc $UniqRootName-AdapterRemovalParameters -f html -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -o $UniqRootName-AdapterRemovalParameters.pdf
        pdfunite $UniqRootName.settings-R.txt-All.pdf $UniqRootName.settings-R.txt-Mate1.pdf $UniqRootName.settings-R.txt-Mate2.pdf $UniqRootName.settings-R.txt-Singletons.pdf $UniqRootName.settings-R.txt-Collapsed.pdf $UniqRootName.settings-R.txt-CollapsedTruncated.pdf $UniqRootName.settings-R.txt-Discarded.pdf $UniqRootName.AdapterRemovalDistributionCharts.pdf
    done
done< UniqRootName

#Archive the Adapter Removal settings files
echo -e "\n\n *** ORGANIZING ADAPTER REMOVAL SETTINGS OUTPUT *** \n\n"
mkdir AdapterRemovalSettings
rsync -acv *.settings AdapterRemovalSettings
rsync -acv *AdapterRemovalParameters.pdf AdapterRemovalSettings
rsync -acv *.AdapterRemovalDistributionCharts.pdf AdapterRemovalSettings
rm *AdapterRemovalParameters.pdf
rm *.AdapterRemovalDistributionCharts.pdf
rm *.settings
rm RInputFiles
rm *Singletons.pdf
rm *Mate1.pdf
rm *Mate2.pdf
rm *Discarded.pdf
rm *Collapsed.pdf
rm *CollapsedTruncated.pdf
rm *All.pdf
rm *R.txt
rm *settings-parameters
rm *-AdapterRemovalParameters

#Run MultiQC
echo -e "\n\n *** RUNNING MULTIQC *** \n\n"
module purge
module load anaconda/colsa
conda activate multiqc-1.10.1
multiqc Stats/ fastqcRawData/ AdapterRemovalSettings/ fastqcTrimmedData/
module purge
module load linuxbrew/colsa

#Archive and compress the raw reads to organize the directory and save disc space
echo -e "\n\n *** ORGANIZING QC REPORT OUTPUT *** \n\n"
rsync -acv *_R1_fastqc.pdf fastqcRawData/
rsync -acv *_R2_fastqc.pdf fastqcRawData/
rsync -acv *.truncated_fastqc.pdf fastqcTrimmedData/
rm *_R1_fastqc.pdf
rm *_R2_fastqc.pdf
rm *.truncated_fastqc.pdf
mkdir QCReports
mkdir QCReports/MultiQCReports
mkdir QCReports/FastQCReports
rsync -acv fastqcRawData QCReports/FastQCReports/
rsync -acv fastqcTrimmedData QCReports/FastQCReports/
rsync -acv multiqc_data QCReports/MultiQCReports
rsync -acv multiqc_report.html QCReports/MultiQCReports
rm -r fastqcRawData
rm -r fastqcTrimmedData
rm -r multiqc_data
rm multiqc_report.html

#Create a list of files to be processed (DOES include Undetermined fastq files)
echo "\n\n *** CALCULATING FILE MD5SUMs *** \n\n"
#Record the md5sums of the files before compression:
#Set date & time stamp variable
current_time=$(date "+%Y.%m.%d-%H.%M")
#Set source directory (i.e., the original directory of files you were backing up)
SourceLocation=$(pwd)
SourceFilename="Premise-$Directory-md5sum-$current_time.txt"
find $SourceLocation -type f \( ! -path '*/.*' \) \( -not -name "$SourceFilename" \)  \( -not -name '*FilesTo*' \)  \( -not -name RInputFiles \)  \( -not -name TrimmedFiles \)  \( -not -name '*RootName' \) -exec md5sum '{}' \; > $SourceFilename

#Archive and compress the raw sequencing reads and reports
echo "\n\n *** ARCHIVING & COMPRESSING RAW SEQUENCING READS & ADAPTER REMOVAL FILES *** \n\n"
ls *.f*q* > AllFilesToArchive
echo Reports >> AllFilesToArchive
echo Stats >> AllFilesToArchive
echo AdapterRemovalSettings >> AllFilesToArchive
tar cvfj "$Directory"_RawSequencingReads.tar.bz2 -T AllFilesToArchive
#find $SourceLocation -type f \( -name $Directory_RawSequencingReads.tar.bz2 \) -exec md5sum '{}' \; >> $SourceFilename
#Backup newly created files
Path=$(pwd)
export https_proxy=http://premise.sr.unh.edu:3128
rclone copy $LocalConfigID:$Path/AdapterRemovalSettings $BoxConfigID:$BoxPath$Directory/AdapterRemovalSettings
$RecycleBin -r Reports
$RecycleBin -r Stats
$RecycleBin AdapterRemovalSettings
$RecycleBin *.f*q*

#Generate LaTeX code for md5sums file
touch $Directory-ELNLaTeXCode-md5sums.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCode-md5sums.txt
echo '\phantomsection' >> $Directory-ELNLaTeXCode-md5sums.txt
echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCode-md5sums.txt
echo '\begin{landscape}' >>  $Directory-ELNLaTeXCode-md5sums.txt
echo "\includepdf[pages=-,frame,scale=0.85,angle=-270,pagecommand={\pagestyle{fancy}}]{$ELNPath$Directory/$SourceFilename.pdf}" >> $Directory-ELNLaTeXCode-md5sums.txt
echo '\end{landscape}' >> $Directory-ELNLaTeXCode-md5sums.txt

#Archive and compress the fastqc data to organize directory and save disc space
echo -e "\n\n *** ARCHIVING & COMPRESSING QC FILES *** \n\n"
rclone copy $LocalConfigID:$Path/QCReports $BoxConfigID:$BoxPath$Directory/QCReports
tar cvfj "$Directory"_QCReports.tar.bz2 QCReports

echo -e "\n\n *** FINAL FILE ORGANIZATION, BACKUP & CLEANUP *** \n\n"
mkdir RawReads-and-SequencingPreProcessingData
rsync -acv *_QCReports.tar.bz2 RawReads-and-SequencingPreProcessingData/
rsync -acv *_RawSequencingReads.tar.bz2 RawReads-and-SequencingPreProcessingData/
$RecycleBin *_QCReports.tar.bz2
$RecycleBin *_RawSequencingReads.tar.bz2

#Add the md5sums of the 2 .tar.bz2 compressed files to the md5sum file:
find $SourceLocation -type f \( -path '*/RawReads-and-SequencingPreProcessingData/*' \)  \( ! -path '*/.*' \) -exec md5sum '{}' \; >> $SourceFilename
#Sort the md5sum file alphabetically by path & filename & add a LaTeX header & footer for formatting purposes:
sed -i 's,'"$SourceLocation"',,g' $SourceFilename
sort -u -k2 $SourceFilename > Formatted-$SourceFilename
while read line; do
    echo $line >> AFormatted-$SourceFilename
    echo ""  >> AFormatted-$SourceFilename
done< Formatted-$SourceFilename
echo -e '---' > Header
echo -e "title: "$SourceFilename"" >> Header
echo -e '---' >> Header
echo -e '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ {.yaml}' >> Header
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" > Footer
cat Header AFormatted-$SourceFilename Footer > Final-$SourceFilename
fold -s -w 135 Final-$SourceFilename > Temp-$SourceFilename.txt
pandoc Temp-$SourceFilename.txt -t latex -V geometry:landscape -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -o $SourceFilename.pdf
rm Temp-$SourceFilename.txt

#Backup additional newly created files
rclone copy $LocalConfigID:$Path --include "*.truncated.bz2" $BoxConfigID:$BoxPath$Directory/QualityTrimmedReads
rclone copy $LocalConfigID:$Path --include "$SourceFilename.pdf" $BoxConfigID:$BoxPath$Directory
rclone copy $LocalConfigID:$Path --include "*ELNLaTeXCode*" $BoxConfigID:$BoxPath$Directory/ELNLaTeXCodeFiles
rsync -acv $SourceFilename* RawReads-and-SequencingPreProcessingData/
mkdir RawReads-and-SequencingPreProcessingData/ELNLaTeXCode
rsync -acv *ELNLaTeXCode* RawReads-and-SequencingPreProcessingData/ELNLaTeXCode

#Cleanup remaining tmp files
rm FilesToProcess
rm RootName
rm UniqRootName
rm AllFilesToProcess
rm AllFilesToArchive
rm TrimmedFiles
rm $SourceFilename*
rm AFormatted-$SourceFilename
rm Formatted-$SourceFilename
rm Final-$SourceFilename
rm Header
rm Footer
$RecycleBin QCReports
rm *ELNLaTeXCode*

#Recommended next step is to start mapping pipeline
#see BWA-MEM_Alignment.sh 
