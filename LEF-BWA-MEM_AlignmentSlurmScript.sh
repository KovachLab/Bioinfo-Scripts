#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --job-name BWA-MEM
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Lindsey.Fenderson@unh.edu

#Script for submitting BWA-MEM_Alignment script to the slurm scheduler, writing all workflow metadata to the slurm output, renaming the slurm output informatively (i.e., with the name of the script run, the directory the script was run on, and the date and time stamp of the run), formatting the file into a nice pdf for input into electronic lab notebook, writing associated LaTeX code to facilitate this, and backing up the slurm output file and pdf to the cloud.

#To reuse this script, confirm the path/script name to be run in lines 27, 30, 35 and 36.

#Set date & time stamp variable
current_time1=$(date "+%Y.%m.%d-%H.%M")
source ${ParamFile}

#Automatically get directory name for the output archive file names
Directory=$(echo $InputDataPath | rev | cut -d'/' -f 2 | rev)
echo -e '---'
echo -e "title: "BWA-MEM Alignment Slurm Output for $Directory""
echo -e '---'
echo -e '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ {.yaml}'
echo -e "The date and run time of this script was:" $current_time1
echo
echo "The script being run, from the following location, and the script's last modified time stamp:"
ls -lh /mnt/lustre/mel/shared/Scripts/BWA-MEM_AlignmentV1.3.sh
echo
echo -e "The header of the script run was: \n"
head /mnt/lustre/mel/shared/Scripts/BWA-MEM_AlignmentV1.3.sh
echo
echo -e "The script was run on the files in the directory: \n"
echo $InputDataPath
echo
echo -e "Running script BWA-MEM_AlignmentV1.3.sh: \n"
/mnt/lustre/mel/shared/Scripts/BWA-MEM_AlignmentV1.3.sh
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
cd ~/BioinformaticRecords
cp ~/BioinformaticScripts/$SLURM_JOB_NAME.$SLURM_JOB_ID.out BWA-MEM_Alignment-$Directory-slurmOutput-$current_time1.yaml
module purge
module load linuxbrew/colsa
sed -i 's/━//g' BWA-MEM_Alignment-$Directory-slurmOutput-$current_time1.yaml

pandoc BWA-MEM_Alignment-$Directory-slurmOutput-$current_time1.yaml -s --highlight-style pygments -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o BWA-MEM_Alignment-$Directory-slurmOutput-$current_time1.pdf
/mnt/lustre/mel/shared/Scripts/rm ~/BioinformaticScripts/$SLURM_JOB_NAME.$SLURM_JOB_ID.out

#Generate LaTeX Code for slurm output file
#Parse directory name & start LaTeX code for ELN
SeqDate=$(echo $Directory | cut -d'_' -f 1 )
SeqPlatform=$(echo $Directory | cut -d'_' -f 2 )
SeqPerson=$(echo $Directory | cut -d'_' -f 3 )
SeqDescription=$(echo $Directory | cut -d'_' -f 4 )
touch $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt
echo '\phantomsection' >> $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt
echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt
echo "\includepdf[pages=-,frame,scale=0.85,pagecommand={\pagestyle{fancy}}]{/Users/Lindsey/Box/KovachLab/Data/Sparrows/GECO/Data/Mapped_Data/$Directory/BioinformaticRecords/BWA-MEMAlignment-$Directory-slurmOutput-$current_time1.pdf}" >> $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt

#Backup newly created files
CurrentPath=$(pwd)
BoxPath='/KovachLab/Data/Sparrows/GECO/Data/Mapped_Data/'
export https_proxy=http://premise.sr.unh.edu:3128
rclone copy Premise:$CurrentPath --include "$Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt" Box:$BoxPath$Directory/ELNLaTeXCodeFiles
rclone copy Premise:$CurrentPath --include "BWA-MEM_Alignment-$Directory-slurmOutput-$current_time1.pdf" Box:$BoxPath$Directory/BioinformaticRecords
/mnt/lustre/mel/shared/Scripts/rm $Directory-ELNLaTeXCode-BWA-MEMAlignmentSlurmOutput.txt
