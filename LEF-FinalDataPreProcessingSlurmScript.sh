#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --job-name SeqPreProcess
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Lindsey.Fenderson@unh.edu

#Script for submitting DataPreProcessingScript to the slurm scheduler, writing all workflow metadata to the slurm output, renaming the slurm output informatively (i.e., with the name of the script run, the directory the script was run on, and the date and time stamp of the run), formatting the file into a nice pdf for input into electronic lab notebook, writing associated LaTeX code to facilitate this, and backing up the slurm output file and pdf to the cloud.

#To reuse this script, be sure to update the working directory of the files that will be processed in line 14, and confirm the path/script name to be run in lines 24, 27, 33 and 34.
source /mnt/lustre/mel/shared/Scripts/DataPreProcessingScriptv1.3.1-Paramfile.txt

module purge
module load linuxbrew/colsa

#Set date & time stamp variable
current_time1=$(date "+%Y.%m.%d-%H.%M")

cd /mnt/lustre/mel/shared/GECO/Data/SparrowRawData/20200728_INVS-SP_LFe_SparrowWholeGenomeShotgun
#Automatically get directory name for the output archive file names
Directory=$(pwd | rev | cut -d'/' -f 1 | rev)
echo -e '---'
echo -e "title: "Data PreProcessing Slurm Output for $Directory""
echo -e '---'
echo -e '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ {.yaml}'
echo -e "The date and run time of this script was:" $current_time1
echo
echo "The script being run, from the following location, and the script's last modified time stamp:"
ls -lh /mnt/lustre/mel/shared/Scripts/DataPreProcessingScriptV1.3.1.sh
echo
echo -e "The header of the script run was: \n"
head /mnt/lustre/mel/shared/Scripts/DataPreProcessingScriptV1.3.1.sh
echo
echo -e "The script was run on the files in the directory: \n"
pwd
echo

echo -e "Running script DataPreProcessingScriptV1.3.1.sh: \n"
/mnt/lustre/mel/shared/Scripts/DataPreProcessingScriptV1.3.1.sh
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
Directory=$(pwd | rev | cut -d'/' -f 1 | rev)
echo $Directory
echo $current_time1
cd ~/BioinformaticRecords
cp ~/BioinformaticScripts/slurm-*.out DataPreProcessing-$Directory-slurmOutput-$current_time1.yaml
module purge
module load linuxbrew/colsa
sed -i 's/━//g' DataPreProcessing-$Directory-slurmOutput-$current_time1.yaml

pandoc DataPreProcessing-$Directory-slurmOutput-$current_time1.yaml --verbose -s --highlight-style pygments -V geometry:"top=2cm, bottom=1.5cm, left=1.5cm, right=1cm" -t latex -o DataPreProcessing-$Directory-slurmOutput-$current_time1.pdf
/mnt/lustre/mel/shared/Scripts/rm ~/BioinformaticScripts/slurm-*.out

#Generate LaTeX Code for slurm output file
#Parse directory name & start LaTeX code for ELN
SeqDate=$(echo $Directory | cut -d'_' -f 1 )
SeqPlatform=$(echo $Directory | cut -d'_' -f 2 )
SeqPerson=$(echo $Directory | cut -d'_' -f 3 )
SeqDescription=$(echo $Directory | cut -d'_' -f 4 )
touch $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt
echo "\subsection{$SeqDate\_$SeqPlatform\_$SeqPerson\_$SeqDescription}" >> $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt
echo '\phantomsection' >> $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt
echo "\addcontentsline{toc}{subsubsection}{$SampleID\_$Species\_$SampleLoc\_$CaptureDate}" >> $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt
echo "\includepdf[pages=-,frame,scale=0.85,pagecommand={\pagestyle{fancy}}]{/Users/Lindsey/Box/KovachLab/Data/Sparrows/GECO/Data/SparrowRawData/$Directory/BioinformaticRecords/DataPreProcessing-$Directory-slurmOutput-$current_time1.pdf}" >> $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt

#Backup newly created files
#RCLONE_CONFIG=/mnt/lustre/mel/leq29/.config/rclone/rclone.conf
#export RCLONE_CONFIG
CurrentPath=$(pwd)
BoxPath='/KovachLab/Data/Sparrows/GECO/Data/SparrowRawData/'
export https_proxy=http://premise.sr.unh.edu:3128
rclone copy Premise:$CurrentPath --include "$Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt" Box:$BoxPath$Directory/ELNLaTeXCodeFiles
rclone copy Premise:$CurrentPath --include "DataPreProcessing-$Directory-slurmOutput-$current_time1.pdf" Box:$BoxPath$Directory/BioinformaticRecords
#rm $Directory-ELNLaTeXCode-DataPreProcessingSlurmOutput.txt
