# A scrip efficiently run and compile code for the hpc.

# To run the scrip you have to change the first 3 lines of code to your own
# needs. The mLine has to be set to the name of the .par file.
# The command variable has to be set to the compiling commmand on the server.
# The NewLocation has to be set to the main location to store the output.
#
# It is important that all the files can be found by the sript, so they have to
# be stored in a ./code subdirectory
#
# After the you get your data it is important to remove them from the server,
# This feature is not yet included in this scrip, and it will not be done
# automatticaly.
#
# Settings for the scrip
mLine='amrvac.par'      # Name of the par file.
command='$AMRVAC_DIR/setup.pl -d=22 -g=16,16 -p=hd -arch=inteldebug'  # compiling code
NewLocation="/scratch/leuven/311/vsc31101/" # Location on the server to store the files

if [ -z "$1" ]
then
echo "No argument supplied"
echo "Usage: "
echo "first argument: Name of the folder created on the server"
echo "Optinal second argument: get and qsub"
echo "    when the extra argument is get the sript will start downloading"
echo "    when the extra argument is qsub th script will resubmit the job"
echo "Example: ./submit.sh rt3d get"
exit 1
fi

echo $1 >> server.log

# Sets the title bar of the terminal
echo -n -e "\033]0;Creating $1\007"

# Creating some aid var
Location=$NewLocation$1
code="code"
S="/"
# Creates the folders
ssh -t -t hpc << EOF
mkdir -p "$Location"
mkdir -p "$Location/run"
mkdir -p "$Location/code"
exit
EOF


if [ -z "$2" ]      # If there is no second argument given, it creates the job
then
#Creating the job file
touch job
cat > job <<EOL
#!/bin/bash -l
#PBS -l walltime=01:00:00
#PBS -l nodes=5:ppn=20
#PBS -M bob.vergauwen@student.kuleuven.be
#PBS -A lp_edu_amrvac_2014
source .bashrc
cd $Location/code
cp ./amrvac ../run/
cp ./$mLine ../run/
cd ../run
mpirun -np 100 ./amrvac -i ./$mLine
EOL

# Copy all the files to the server
sftp  hpc << EOF
put -r code/ $Location
put job $Location
exit
EOF
# Compling and running submiting the code
ssh -t -t hpc << EOF
cd  "$Location/code"
$command
make
cd ../
qsub -N "$1" job
exit
EOF
fi

if [ "$2" == "qsub" ] # If the second argument is qsub, just resubmit the code
then
# submitting the file
ssh -t -t hpc << EOF
cd  "$Location"
qsub -N "$1" job
exit
EOF
fi

if [ "$2" == "get" ] # If the second argument is get, compress the files and download them
then
# Getting the files
echo 'Getting the file'
ssh -t -t hpc << EOF
cd  "$Location"
zip -r ./data.zip ./run/*.vtu
exit
EOF
mkdir -p data
sftp  hpc << EOF
get $Location/data.zip ./data/data.zip
exit
EOF
fi
