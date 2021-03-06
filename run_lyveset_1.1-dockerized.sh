#! /bin/bash
#Author: Logan Fink
#Usage: script to select a reference genome for input to lyve-SET and run lyve-SET
#Permission to copy and modify is granted without warranty of any kind
#Last revised 03/12/18

run_quit=false
while getopts :q opt; do
    case $opt in
        q)run_quit=true ;;
    esac
done

shift "$(( OPTIND-1 ))"

# set number of max number of threads available
THREADS=$(nproc --all)
echo "Number of threads set to: $THREADS"
export THREADS

##### Check to see if docker is installed #####
if [ -z $(which docker) ]; then
   echo 'Docker is not installed, Please see https://github.com/StaPH-B/scripts/blob/master/image-information.md#docker-ce for instrucions on how to install Docker. Now exiting.'
   exit
else
    echo "$(docker --version) is installed."
fi

##### Function and check to see if Docker images are downloaded, if not, download them with docker pull #####
docker_image_check () {
if [ -z $(docker images -q $1) ]; then
    docker pull $1
else
    echo "Docker image $1 already exists locally."
fi
}
# gotta check em all! gotta check em all! Dock-er-mon!
echo "Now checking to see if all necessary docker images are downloaded..."
docker_image_check staphb/lyveset:2.0.1


##### Create the snp_counts directory #####
python /home/staphb/scripts/pairwise_matrix_dist*.py ./roary/core_gene_alignment.aln
ref=$(head -1 ./snp_counts/avg_snp_differences_ordered | cut -d ' ' -f 1)

##### Create lyveset directory. If it exists, delete it and create it #####
rm -rf lyveset 2>/dev/null
docker run --rm=True -u $(id -u):$(id -g) -v $PWD:/data staphb/lyveset:2.0.1 \
set_manage.pl --create lyveset

cp ./spades_assembly_trim/*$ref*/contigs.fasta lyveset/reference/$ref.fasta
cp ./clean/*.fastq.gz lyveset/reads/
docker run --rm=True -u $(id -u):$(id -g) -v $PWD:/data staphb/lyveset:2.0.1 \
launch_set.pl lyveset -ref lyveset/reference/$ref.fasta --allowedFlanking 5 --min_alt_frac 0.95 --min_coverage 20 --numcpus ${THREADS}

if $run_quit; then
    sudo shutdown -h now
fi
