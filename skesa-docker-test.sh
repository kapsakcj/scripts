#!/bin/bash

#This function will check if the file exists before trying to remove it
remove_file() {
    if [ -e $1 ];then
        echo "File $1 has been removed"
        rm -rf $1
    fi
}

#This function will check to make sure the directory doesn't already exist before trying to create it
make_directory() {
    if [ -e $1 ]; then
        echo "Directory "$1" already exists"
    else
        mkdir $1
        echo "Directory "$1" has been created"
    fi
}

if [ -z "$SEQUENCE_LEN" ]; then
    SEQUENCE_LEN=5000000
fi
echo "Sequence Length: $SEQUENCE_LEN"

#####Number of threads-specific to linux
THREADS=$(nproc --all)
echo "Number of threads set to: $THREADS"
export THREADS

##### Move all fastq files from fastq_files directory up one directory, remove fastq_files folder #####
if [[ -e ./fastq_files ]]; then
    echo "Moving fastq files from ./fastq_files to ./ (top-level DIR)"
    mv ./fastq_files/* .
    rm -rf ./fastq_files
fi

declare -a srr=() #PASTE IN ANY SRR NUMBERS INTO FILE named: SRR
while IFS= read -r line; do
    srr+=("$line")
done < ./SRR 2>/dev/null
#find . -maxdepth 1 -name '*fastq*' |cut -d '-' -f 1|cut -d '_' -f 1 |cut -d '/' -f 2 >tmp1 #output all fastq file identifiers in cwd to file tmp1 (the delimiters here are '-' and '_')
find . -maxdepth 1 -name '*fastq*' |cut -d '_' -f 1 |cut -d '/' -f 2 >tmp1 #output all fastq file identifiers in cwd to file tmp1 (the delimiters here are '_')
declare -a tmp=()
tmp1='tmp1'
tmpfile=`cat $tmp1`
for line in $tmpfile; do
    tmp=("${tmp[@]}" "$line");
done
id=("${tmp[@]}" "${srr[@]}") #Combine the automatically generated list with the SRR list
id=($(printf "%s\n" "${id[@]}"|sort -u)) #Get all unique values and output them to an array
echo ${id[@]}
remove_file tmp

##### Fetch and fastq-dump all reads from NCBI identified by "SRR" #####
for i in ${srr[@]}; do
    echo $i
    if [[ -n "$(find *$i* 2>/dev/null)" ]]; then
        echo "Files are here."
    else
        export i
        echo 'prefetching '$i'...'
        print_next_command $LINENO ${i}
        docker run -e i --rm=True -u $(id -u):$(id -g) -v $PWD:/data staphb/sratoolkit:2.9.2 /bin/bash -c \
        'prefetch -O /data '${i}
        echo 'now running fasterq-dump in container'
        print_next_command $LINENO ${i}
        docker run -e THREADS -e i --rm=True -u $(id -u):$(id -g) -v $PWD:/data staphb/sratoolkit:2.9.2 /bin/bash -c \
        'fasterq-dump --skip-technical --split-files -t /data/tmp-dir -e ${THREADS} -p '${i}'.sra'
        mv ${i}.sra_1.fastq ${i}_1.fastq
        mv ${i}.sra_2.fastq ${i}_2.fastq
        pigz ${i}_1.fastq
        pigz ${i}_2.fastq
        rm ${i}.sra
    fi
done
remove_file tmp-dir

##### These are the QC trimming scripts as input to trimClean #####
make_directory clean
echo "cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning cleaning"
for i in *R1_001.fastq.gz; do
    b=`basename ${i} _R1_001.fastq.gz`;
    if [[ -n "$(find -path ./clean/${b}.cleaned.fastq.gz)" ]]; then
        continue
    else
        echo "LYVESET CONTAINER RUNNING SHUFFLEREADS.PL"
        print_next_command $LINENO ${b}
        docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/lyveset:2.0.1 \
        run_assembly_shuffleReads.pl /data/${b}"_R1_001.fastq.gz" /data/${b}"_R2_001.fastq.gz" > clean/${b}.fastq;
        echo ${b};
        echo "LYVESET CONTAINER RUNNING TRIMCLEAN.PL"
        print_next_command $LINENO ${b}
        docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/lyveset:2.0.1 \
        run_assembly_trimClean.pl --numcpus ${THREADS} -o /data/clean/${b}.cleaned.fastq.gz -i /data/clean/${b}.fastq --nosingletons;
        remove_file clean/${b}.fastq;
    fi
done
if [ -s "SRR" ]; then
    for j in *_1.fastq.gz; do
        c=`basename ${j} _1.fastq.gz`;
        if [[ -n "$(find -path ./clean/${c}.cleaned.fastq.gz)" ]]; then
            continue
        else
            echo "(run_assembly_shuffleReads.pl)Interleaving reads for:"${c}" using lyveset docker container"
            print_next_command $LINENO ${c}
            docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/lyveset:2.0.1 \
            run_assembly_shuffleReads.pl /data/${c}"_1.fastq.gz" /data/${c}"_2.fastq.gz" > clean/${c}.fastq;
            echo "(run_assembly_trimClean.pl) Trimming/cleaning reads for:"${c}" using lyveset docker container";
            print_next_command $LINENO ${c}
            docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/lyveset:2.0.1 \
            run_assembly_trimClean.pl --numcpus ${THREADS} -i /data/clean/${c}.fastq -o /data/clean/${c}.cleaned.fastq.gz --nosingletons;
            remove_file clean/${c}.fastq;
        fi
    done
else
    echo "There are no SRR numbers in this run"
fi
rm ./clean/\** 2>/dev/null

##### Run SKESA de novo genome assembler on all cleaned, trimmed, fastq files #####
make_directory ./skesa
for i in ${id[@]}; do
    make_directory skesa/${i}
    if [[ -n "$(find -path ./skesa/$i/contigs.fasta 2>/dev/null)" ]]; then #This will print out the size of the spades assembly if it already exists
        size=$(du -hs ./skesa/$i/$i.skesa.contigs.fasta | awk '{print $1}');
        echo 'File exists and is '$size' big.'
    else
        echo 'constructing skesa assemblies for '$i', could take some time...'
        # exporting `i` variable to make it available to the docker container
        export i
        docker run -e i -e THREADS --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/skesa:2.3.0 /bin/bash -c \
        'skesa --fastq /data/clean/*$i*.cleaned.fastq.gz --use_paired_ends --cores ${THREADS} --contigs_out /data/skesa/${i}/${i}.skesa.contigs.fasta'
    fi
done

##### Run quast assembly statistics for the skesa assemblies
make_directory quast-skesa
for i in ${id[@]}; do
     if [[ -n "$(find -path ./quast-skesa/${i}_output_file 2>/dev/null)" ]]; then
        echo "Skipping ${i} - It's skesa assembly has already been QUASTed."
    else
    	docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/quast:5.0.0 \
    	quast.py --fast -t ${THREADS} /data/skesa/$i/$i.skesa.contigs.fasta -o /data/quast-skesa/$i
    fi
done

#### Remove the tmp1 file that lingers #####
remove_file tmp1

##### Move all of the fastq.gz files into a folder #####
make_directory fastq_files
for i in ${id[@]}; do
    if [[ -n "$(find *$i*fastq.gz)" ]]; then
        mv *$i*fastq.gz fastq_files
    fi
done
todays_date=$(date)
echo "*******************************************************************"
echo "skesa-docker-test pipeline has finished on "$todays_date"."
echo "*******************************************************************"
