#!/usr/bin/env python

# ./wrapper.py -i SRR821573.cram -r 0.1:0.9:0.1 -a /scratch0/genomes/hg38/annotation/hg38_p8.biotype_flt.cls.gff3 -e /scratch0/genomes/hg38/hg38.fa

import os
import argparse
import sys
import subprocess

# This function shall calculate the correct number of:
# 1. forks to maintain simultaneously
# 2. threads to allocate to each fork
def CalcNumThreads(requestedTreads):
    cpuinfo = open('/proc/cpuinfo','r')
    forkNum = 0
    for line in cpuinfo:
        try:
            if line.split()[0] == "siblings":
                forkNum = line.split()[2]/2/requestedTreads
                if forkNum == 0:
                    forkNum = 1
        except:
            pass
    cpuinfo.close()
    return forkNum # needs to be something else

def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1

# This function will be used in order to verify whether the input information is submitted in bulk or as a single piece and process the information accordingly
def checkDirFile(inStr):
    toAnalyze = []
    rootDir = os.popen("pwd").read()[:-1]

    if os.path.isdir(inStr):
        for d in os.listdir(inStr):
            if not os.path.isdir(d):
                # is directory. compile a list of all files with full paths
                toAnalyze.append(os.path.abspath(inStr+"/"+d))
            else:
                # Needs to raise an exception
                return "ERROR"
        return toAnalyze
    else:
        return [inStr]

def main(argv):

    curPath = os.path.dirname(os.path.realpath(__file__))
    outDir = curPath

    threads = 1

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
    parser.add_argument('-f','--forks',type=int,help="Indicate the maximum numberof forks to maintain at the same time. Each fork will receive either precalculated or predefined by -t number of threads")
    parser.add_argument('-o','--out',type=str,help="Directory where all the output will be saved")
    # possible need to specify output path

    args=parser.parse_args()
    inputs = args.input
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

    forkNum = CalcNumThreads(threads)

    if args.stats is not None:
        open(outLog,'a').close()
        print("Your stats log is created in"+outLog)
    else:
        open(outDir+"/stats.log",'a').close()
        print("The default log was created in"+outDir+"/stats.log")

    if args.out is not None and os.path.isdir(args.out):
        outDir = args.out

    # For MARCC need to add parallelization within this wrapper:
    # this level of paralellization will allow running the wrapper simultaneously on multiple tissue samples (read as input files).
    # perhaps do it simply by system forking the process - is very likely to be the most efficient way of doing it

    def child(path):
    # Child should receive its own pseudo-terminal for proper threading of the application
        # Add parsing of the filename
        baseDirName = path.split("/")[-1].split(".")[:-1][0]

        if(not os.path.exists(outDir+"/downsamp/"+baseDirName)):
            os.makedirs(outDir+"/downsamp/"+baseDirName)
        for i in xfrange(covRange[0],covRange[1],covRange[2]):
            if(os.path.exists(outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")):
                os.system("rm -r "+outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")

        if(not os.path.exists(outDir+"/assembly/"+baseDirName)):
            os.makedirs(outDir+"/assembly/"+baseDirName)
        for i in xfrange(covRange[0],covRange[1],covRange[2]):
            if(os.path.exists(outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))):
                os.system("rm -r "+outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))
            os.makedirs(outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))

        if(not os.path.exists(outDir+"/statsAl/"+baseDirName)):
            os.makedirs(outDir+"/statsAl/"+baseDirName)
        for i in xfrange(covRange[0],covRange[1],covRange[2]):
            if(os.path.exists(outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))):
                os.system("rm -r "+outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))
            os.makedirs(outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))

        # Perhaps also add a step to break the alignment into functional pieces based on the reference, so that separate directories/files may be created for different genes. Analysis will be much faster and perhaps easier

        # also do not forget to create a separate directory for a base assembly
        # perhaps just save it in the root of each directory where all the subdirs for downsampled assemblies are

        #==============================================
        #first run the base assembly as described above
        #==============================================


        # print("==================================")
        # print(sequenceRef)
        # print("==================================")
        assembleFullCov = "samtools view -h -T "+os.path.abspath(sequenceRef)+" "+path+" | stringtie -p "+str(threads)+" -m 150 -G "+os.path.abspath(annotationRef)+" -o "+outDir+"/assembly/"+baseDirName+"/"+baseDirName+".gtf -A "+outDir+"/assembly/"+baseDirName+"/Genes.tab -"
        os.system(assembleFullCov)

        # Now lets write the first downsampling method using samtools
        for i in xfrange(covRange[0],covRange[1],covRange[2]):
            downSampleCmd = "samtools view -h -s "+str(i)+" -C "+path+" > "+outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram"
            os.system(downSampleCmd)


        # for each child Fork create a temporary log in its location so that information can be saved sequentially
        # At the end of all the child processes unify all the temp logs into a single comprehensive under the parent thread

        # Exit after the process is over
        os._exit(0)

    def parent(inputs):
        # the block below is responsible for creating the base directories
        # Again, need to integrate an output parameter passing and if specified create basedir in accordance with the output parameter
        if(not os.path.exists(outDir+"/downsamp")):
            os.makedirs(outDir+"/downsamp/")
        if(not os.path.exists(outDir+"/assembly")):
            os.makedirs(outDir+"/assembly/")
        if(not os.path.exists(outDir+"/statsAl")):
            os.makedirs(outDir+"/statsAl/")

        theorizedForkNum = CalcNumThreads(threads) # Calculates the number of forks to make
        if args.forks is not None and args.forks <= theorizedForkNum*2: # safeguards against excessive forking
            forkNum = args.forks
        else:
            forkNum = theorizedForkNum

        inputAls = checkDirFile(inputs)
        childPIDS = [] # list to be populated with the fork PIDS
        # while ( we parse through the list of tissue alignments to analyze )
        while(len(inputAls) != 0): # need to creating waiting for the process PIDs to finish before forking again
            newpid = os.fork()
            path = inputAls.pop()
            if newpid == 0:
                child(os.path.abspath(path))
            else:
                childPIDS.append((os.getpid(),newpid))
                print("PARENT: NEW CHILD CREATED: ",childPIDS)
                if len(childPIDS) >= forkNum:
                    expiredPID = os.wait()
                    print("EXPIRED: ",expiredPID)
                    childPIDS.remove((os.getpid(),expiredPID[0]))
                    
            while(len(childPIDS) > 0):
                expiredPID = os.wait()
                print("FINAL EXPIRATION: ",expiredPID)
                childPIDS.remove((os.getpid(),expiredPID[0]))

        # The parent currently exits before the child ends the process. Parent needs to wait for all childrent to finish

    parent(inputs)

# also important to create a safeguard for threading:
# if the user specifies 8 threads to be used, the script meeds to analyze the max available and how muchforking will be done and calculate the number of threads to be allocated per fork (read for stringtie)

# When producing plots the following idea might come very useful:
# 1. Why not use K-means clustering in order to cluster genes by similarity in terms of TPM/coverage/"other parameters". Otherwise selecting for genes or transcripts might not be possible.
# With K-means it will also be possible to bring some insight into what is changing with downsampling and what is being affected

if __name__=="__main__":
    main(sys.argv[1:])
