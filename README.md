# type_pipe-dockerized

### Requirements
  * `docker`
  * `pigz`
  * `~/github-repos/scripts` must be in PATH. Clone git repo here.
  

### Usage
Requires either :
  * a text file named: `SRR` containing SRR ID's in a list:
```
SRR1234567
SRR1234568
[...]
```
  * PE fastq files that end in either `R1_001.fastq.gz` or `_1.fastq.gz`
  
To run:
```bash
cd dir-containing-fastq-or-SRR-file
run_type_pipe_2.4-dockerized.sh
```

This pipeline takes in raw fastq files, cleans them, runs them through kraken, mash, SPAdes, QUAST, SeqSero1 & 2, SISTR, Serotypefinder, and abricate for antibiotic resistance determination.



