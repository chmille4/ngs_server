#Ultralightweight NGS Server

   Converts BAM and VCF files to json for consumption by a web application

 
## install for mac and linux (no windows)
1. download and install cmake

2. run "gem install ngs_server" in the terminal

## usage
      # in terminal
      ngs_server --help
      
      # add data directories to be served
      ngs_server add path-to-dir
      
      # start server
      ngs_server start
      ngs_server start --port 3000
      # daemonize
      ngs_server start -d
      
      # stop server if daemonized
      ngs_server stop
      

## access data
      
see available data sources

* http://0.0.0.0:4569/sources/bam
* http://0.0.0.0:4569/sources/vcf
      
<br/>invoke file: hg18.bam with coordinates 1 to 100000 on chromosome 22

* http://0.0.0.0:4569/json/bam/hg.18?min=1&max=10000&segment=22
      
<br/>vcf files must be compressed with index (see tabix)<br/>
invoke file: genotypes.vcf.gz with coordinates 1073361 to 1238825 on chromosome 1

* http://0.0.0.0:4569/json/vcf/genotypes.vcf.gz?segment=1&min=1073361&max=1238825


## Example
    # create data directory
    mkdir example-data && cd example-data
    
    # download bam file and index
    wget http://bioinformatics.bc.edu/marthlab/download/small.bam
    wget http://bioinformatics.bc.edu/marthlab/download/small.bam.bai
    
    # add data directory to ngs_server to be served
    cd ../
    ngs_server add example-data/
    
    # start server
    ngs_server start -p 4569
    
    # paste following url in browser
    # to download region of bam as json in web browser
    # http://0.0.0.0:4569/json/bam/example-data/small.bam?min=14596532&max=14699000&segment=22
    