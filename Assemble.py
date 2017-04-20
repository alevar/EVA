import os
import argparse
import sys
import random
import multiprocessing
import signal
import glob
import shutil

# The following function is meant to parse the data and produce an output for plotting
def parseDir(path,sample,params):
    assemPath=path+"/transcripts.gtf"
    cmpCMD = """ " """+sample+""" " "\\t" $1 "\\t" $11 "\\t" $10 "\\t" $4 "\\t" $5 "\\t" """+str(params[2])+""" "\\t" """+str(params[1])+""" "\\t" $12 "\\t" $14"""
    awkCMD = """awk -F "\\t" '{print """+cmpCMD+"}' "+assemPath+" >> "+params[0]
    os.system(awkCMD)

# This function shall verify and readjust the correct number of:
# 1. forks to maintain simultaneously
# 2. threads to allocate to each fork
def verifyThreading(requestedForks,requestedThreads):
    if (requestedThreads*requestedForks <= multiprocessing.cpu_count()) and (requestedThreads<=multiprocessing.cpu_count()):
        return (requestedForks,requestedThreads)
    else:
        return(1,1)

def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1

# This function will be used in order to verify whether the input information is submitted in bulk or as a single piece and process the information accordingly
def checkDirFile(inStr):
    toAnalyze = []
    if not isinstance(inStr,basestring):
        return inStr
    if os.path.isdir(inStr):
        for d in os.listdir(inStr):
            if not os.path.isdir(d):
                # is directory. compile a list of all files with full paths
                toAnalyze.append(os.path.abspath(inStr+"/"+d))
            else:
                # Needs to raise an exception
                return "ERROR: check the spelling of your input files and directories"
        return toAnalyze
    else:
        return [inStr]

def child(path,outDir,covRange,numReps,sequenceRef,annotationRef,threads,cont):
    baseDirName = path.split("/")[-1].split(".")[:-1][0]
    finDir = outDir+"/"+baseDirName+"/"+baseDirName
    if not cont==False:
        # If this method with the script works - consider writing the config file from here rather than passing parameters
        for scalefactor in xfrange(cont,covRange[1],covRange[2]):
            randSeed = random.sample(range(1,10000),numReps) # Need to think where to place this line
            for rep in range(numReps):
                scriptCMD = "./rnaseq_al_pipe.sh "+path+" "+finDir+" "+str(randSeed[rep]+scalefactor)+" "+str(rep)+" "+sequenceRef+" "+annotationRef+" "+str(threads)
                os.system(scriptCMD)
    else:
        # If this method with the script works - consider writing the config file from here rather than passing parameters
        for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
            randSeed = random.sample(range(1,10000),numReps) # Need to think where to place this line
            for rep in range(numReps):
                scriptCMD = "./rnaseq_al_pipe.sh "+path+" "+finDir+" "+str(randSeed[rep]+scalefactor)+" "+str(rep)+" "+sequenceRef+" "+annotationRef+" "+str(threads)
                os.system(scriptCMD)
    os.remove(finDir+".sam")
    os._exit(0)

childPIDS = []
def parent(inputs,outDir,covRange,numReps,sequenceRef,annotationRef,tf,cont):
    if not cont==False:
        print("=============Continue=============")
        for sf in xfrange(cont,covRange[1],covRange[2]):
            for directory in glob.glob(outDir+"/*/*"+str(sf)+"*/"):
                shutil.rmtree(directory)

    global childPIDS
    # while ( we parse through the list of tissue alignments to analyze )
    for inputFile in inputs:
        path = os.path.abspath(inputFile)
        if len(childPIDS) >= tf[0]:
                childPIDS[0].join()
                childPIDS.remove(childPIDS[0])
        else:
            p = multiprocessing.Process(target=child, args=(os.path.abspath(path),outDir,covRange,numReps,sequenceRef,annotationRef,tf[1],cont,))
            childPIDS.append(p)
            p.start()

    while(len(childPIDS) > 0):
        childPIDS[-1].join()
        childPIDS.remove(childPIDS[-1])

def main(args):
    inputs = checkDirFile(args.input)
    outDir = "./"
    if args.out is not None:
        outDir = os.path.abspath(args.out)
        os.system("mkdir -p "+outDir)
    else:
        curPath = os.path.dirname(os.path.realpath(__file__))
        outDir = curPath

    assert len(args.range.split(":"))==3
    covRange = [float(args.range.split(":")[0]),float(args.range.split(":")[1]),float(args.range.split(":")[2])]

    numReps = args.random
    sequenceRef = args.reference
    annotationRef = args.annotation
    assert os.path.isfile(annotationRef)
    assert os.path.isfile(sequenceRef)

    threading = verifyThreading(args.forks,args.threads)

    if args.stats is not None:
        statFilePath = args.stats
    else:
        statFilePath = outDir+"/stats.log"

    parent(inputs,outDir,covRange,numReps,sequenceRef,annotationRef,threading,args.cont)

    # here we begin parsing the assembled output .gtf
    statFile=open(statFilePath,'w+')
    statFile.write("@TISSUE CHR REFID TNAME START END SAMPLE FACTOR COV TPM\n")
    statFile.close()

    # [statPath,scalingFactor,repetition,baseCov,baseTPM]
    for inputFile in inputs:
        path = os.path.abspath(inputFile)
        sample = path.split("/")[-1].split(".")[:-1][0]
        for sf in xfrange(covRange[0],covRange[1],covRange[2]):
            for rep in range(numReps):
                params=[statFilePath,sf,rep]
                endDir = outDir+"/"+sample+"/"+sample+"_F:"+str(sf)+"_R:"+str(rep)
                parseDir(endDir,sample,params)