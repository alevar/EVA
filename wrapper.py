#!/usr/bin/env python

# ./wrapper.py -i SRR821573.cram -r 0.1:0.9:0.1 -a /scratch0/genomes/hg38/annotation/hg38_p8.biotype_flt.cls.gff3 -e /scratch0/genomes/hg38/hg38.fa

import os
import argparse
import sys
import subprocess

def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1

def main(argv):

    # For MARCC need to add parallelization within this wrapper:
    # this level of paralellization will allow running the wrapper simultaneously on multiple tissue samples (read as input files).
    # perhaps do it simply by system forking the process - is very likely to be the most efficient way of doing it

    def child():
    # Child should receive its own pseudo-terminal for proper threading of the application
        return

    def parent(inputs):
        # while ( we parse through the list of tissue alignments to analyze )
        while(len(inputs) != 0):
            newpid = os.fork()
            if newpid == 0:
                child()
            else:
                print("hola", os.getpid(),newpid)

    curPath = os.path.dirname(os.path.realpath(__file__))

# Arguments to consider adding in future:
# 1. Directory with all input alignments
# 2. Directory with all input reference annotations
# 3. Directory with all input reference sequences for samtools to be able to view the CRAM files correctly

    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument('-i','--input',required=True,type=str,help="path to the CRAM allignment which is to be used in coverage efficacy analysis by downsampling")
    parser.add_argument('-s','--stats',type=str,help="Path to the location of the informational log about the simulation")
    parser.add_argument('-r','--range',required=True,type=str,help="Specify downsampling range and step where 1.0 is full and 0.0 is none. Format: [start]:[stop]:[step]. Stop point is not included")
    parser.add_argument('-a','--annotation',required=True,type=str,help="Provide the path to the reference annotation")
    parser.add_argument('-e','--reference',required=True,type=str,help="Provide path to the referernce sequence in fasta format which will be used in order to transalate CRAM into SAM and feed it into the stringtie for assembly")
    parser.add_argument('-t','--threads',type=int,help="Indicate the maximum number of threads to be used by the pipeline")

    args=parser.parse_args()
    inCRAM = args.input
    assert len(args.range.split(":"))==3
    covRange = [float(args.range.split(":")[0]),float(args.range.split(":")[1]),float(args.range.split(":")[2])]

    sequenceRef = args.reference
    annotationRef = args.annotation
    assert os.path.isfile(annotationRef)
    assert os.path.isfile(sequenceRef)

    if args.stats is not  None:
        outLog = args.stats

    if args.threads is not None:
        threads = args.threads
    else:
        threads = 1

# Add parsing of the filename
    baseDirName = inCRAM.split("/")[-1].split(".cram")[0]

    # Perhaps also add a step to break the alignment into functional pieces based on the reference, so that separate directories/files may be created for different genes. Analysis will be much faster and perhaps easier

    if(not os.path.exists(curPath+"/downsamp")):
        os.makedirs(curPath+"/downsamp/")
        os.makedirs(curPath+"/downsamp/"+baseDirName)
    for i in xfrange(covRange[0],covRange[1],covRange[2]):
        if(os.path.exists(curPath+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")):
            os.system("rm -r downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")

    if(not os.path.exists(curPath+"/assembly")):
        os.makedirs(curPath+"/assembly/")
        os.makedirs(curPath+"/assembly/"+baseDirName)
    for i in xfrange(covRange[0],covRange[1],covRange[2]):
        if(os.path.exists(curPath+"/assembly/"+baseDirName+"/"+baseDirName+str(i))):
            os.system("rm -r assembly/"+baseDirName+"/"+baseDirName+str(i))
        os.makedirs(curPath+"/assembly/"+baseDirName+"/"+baseDirName+str(i))

    if(not os.path.exists(curPath+"/statsAl")):
        os.makedirs(curPath+"/statsAl/")
        os.makedirs(curPath+"/statsAl/"+baseDirName)
    for i in xfrange(covRange[0],covRange[1],covRange[2]):
        if(os.path.exists(curPath+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))):
            os.system("rm -r statsAl/"+baseDirName+"/"+baseDirName+str(i))
        os.makedirs(curPath+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))

    if args.stats is not None:
        open(outLog,'a').close()
        print("Your stats log is created in"+outLog)
    else:
        open(curPath+"stats.log",'a').close()
        print("The default log was created in"+curPath+"stats.log")

# also do not forget to create a separate directory for a base assembly
# perhaps just save it in the root of each directory where all the subdirs for downsampled assemblies are

#==============================================
#first run the base assembly as described above
#==============================================
    #assembleFullCov = "samtools view -h -T "+sequenceRef" "+inCRAM+" | stringtie -p "+str(threads)+" -m 150 -G "+annotationRef+" -o ./assembly/"+baseDirName+"/"+baseDirName+".gtf -A ./assembly/"+baseDirName+"/Genes.tab -"
    #os.system(assembleFullCov)

    # Now lets write the first downsampling method using samtools
    for i in xfrange(covRange[0],covRange[1],covRange[2]):
        downSampleCmd = "samtools view -h -s "+str(i)+" -C "+inCRAM+" > ./downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram"
        os.system(downSampleCmd)

# also important to create a safeguard for threading:
# if the user specifies 8 threads to be used, the script meeds to analyze the max available and how muchforking will be done and calculate the number of threads to be allocated per fork (read for stringtie)

# When producing plots the following idea might come very useful:
# 1. Why not use K-means clustering in order to cluster genes by similarity in terms of TPM/coverage/"other parameters". Otherwise selecting for genes or transcripts might not be possible.
# With K-means it will also be possible to bring some insight into what is changing with downsampling and what is being affected

if __name__=="__main__":
    main(sys.argv[1:])
