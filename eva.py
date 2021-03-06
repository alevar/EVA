#!/usr/bin/env python

# EVA - Expression Variation Analyzer
# A tool for assessing the amount of variation in tpm estimation as a function of coverage and number of aligned reads

import os
import sys
import argparse
import signal
import SelectCrams
import Assemble
import Analyze
import Plot

def eva(argv):

    parser = argparse.ArgumentParser(description='''Help Page''')
    subparsers = parser.add_subparsers(help='sub-command help')

#===========================================
#==================SELECT===================
#===========================================

    parser_select = subparsers.add_parser('select',
                                help='select help')
    parser_select.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                help="directory which contains log output files from bowtie2")
    parser_select.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./alignmentLog.txt",
                                help="output for the csv containing information about the location of selected alignments")
    parser_select.add_argument('-n',
                                '--number',
                                required=False,
                                type=int,
                                default=20,
                                help="the number of alignments to select for successive analysis")
    parser_select.add_argument('-r',
                                '--range',
                                required=False,
                                type=str,
                                default="auto",
                                help="Range of the number of aligned reads to evaluate. Default auto will optimize for both the greatest mean number of aligned reads and the samllest range. Custom range may be specified as: minBound:maxBound. If a single integer is provided the algorithm will bin data into this many bins and select alignments based on the number of aligned pairs/reads and the requested number of desired alignments.")
    parser_select.set_defaults(func=SelectCrams.main)

#===========================================
#==================ASSEMBLE=================
#===========================================

    parser_assemble = subparsers.add_parser('assemble',
                                help='assemble help')
    parser_assemble.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                nargs="*",
                                help="path to the CRAM allignment which is to be used in coverage efficacy analysis by downsampling")
    parser_assemble.add_argument('-s',
                                '--stats',
                                type=str,
                                help="Path to the location of the informational log about the simulation")
    parser_assemble.add_argument('-r',
                                '--range',
                                required=True,
                                type=str,
                                help="Specify downsampling range and step where 1.0 is full and 0.0 is none. Format: [start]:[stop]:[step]. Stop point is not included")
    parser_assemble.add_argument('-a',
                                '--annotation',
                                required=True,
                                type=str,
                                help="Provide the path to the reference annotation")
    parser_assemble.add_argument('-e',
                                '--reference',
                                required=True,
                                type=str,
                                help="Provide path to the referernce sequence in fasta format which will be used in order to transalate CRAM into SAM and feed it into the stringtie for assembly")
    parser_assemble.add_argument('-t',
                                '--threads',
                                type=int,
                                default=1,
                                help="Indicate the maximum number of threads to be used by the pipeline")
    parser_assemble.add_argument('-f',
                                '--forks',
                                type=int,
                                default=1,
                                help="Indicate the maximum numberof forks to maintain at the same time. Each fork will receive either precalculated or predefined by -t number of threads")
    parser_assemble.add_argument('-o',
                                '--out',
                                type=str,
                                help="Directory where all the output will be saved")
    parser_assemble.add_argument('-w',
                                '--random',
                                required=True,
                                type=int,
                                help="The number of random sampling repetitions to perform for each genomic region at each scaling factor")
    parser_assemble.add_argument('-c',
                                '--cont',
                                required=False,
                                type=float,
                                default=False,
                                help="The value specify at which downsampling factor to resume the run. If set to true - will attempt to resume the previous run (happened a few times that i accidentally exited tmux and killed the process before it was over). It requires that the rest of the parameters are set the same way as before.")
    parser_assemble.add_argument('-l',
                                '--gene',
                                required=False,
                                action="store_true",
                                help="If set true will also output data for gene-level analysis")
    parser_assemble.set_defaults(func=Assemble.main)

#===========================================
#==================ANALYZE==================
#===========================================

    parser_analyze = subparsers.add_parser('analyze',
                                help='analyze help')
    parser_analyze.add_argument('-i',
                                '--input',
                                required=False,
                                type=str,
                                help="path to the stats.log generated by or in the format of eva assemble")
    parser_analyze.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./analysis",
                                help="output directory")
    parser_analyze.add_argument('-c',
                                '--coverage',
                                required=False,
                                type=str,
                                default="full",
                                help="specify the coverage range to analyze and output")
    parser_analyze.add_argument('-t',
                                '--top',
                                required=False,
                                type=int,
                                help="Select only top N transcripts by expression level (TPM) for analysis")
    parser_analyze.add_argument('-d',
                                '--de',
                                required=False,
                                action="store_true",
                                help="conduct estimation of statistically significant TPM values for differential expression")
    parser_analyze.add_argument('-l',
                                '--gene',
                                required=False,
                                type=str,
                                help="provide gene-lvel assembly results for analysis")
    parser_analyze.add_argument('-f',
                                '--fast',
                                action="store_true",
                                help="speed up the grouping and statistical analysis by using grouped by ID as input to sf")
    parser_analyze.set_defaults(func=Analyze.main)

#===========================================
#===================PLOT====================
#===========================================

    parser_plot = subparsers.add_parser('plot',
                                help='plot help')
    parser_plot.add_argument('-r',
                                '--resolution',
                                required=False,
                                type=str,
                                default="15:15",
                                help="resolution of the plots provided as: width:height")
    parser_plot.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./analysis",
                                help="output directory")
    parser_plot.add_argument('-c',
                                '--coverage',
                                required=False,
                                type=str,
                                default="full",
                                help="what range of coverage to plot")
    parser_plot.add_argument('-f',
                                '--sf',
                                type=str,
                                help="plot data grouped by sampling factor")
    parser_plot.add_argument('-i',
                                '--id',
                                type=str,
                                help="plot data grouped by unique ID")
    parser_plot.add_argument('-d',
                                '--de',
                                type=str,
                                help="plot student t test results for diff expression/fold expression increase signifficance")
    parser_plot.add_argument('-g',
                                '--gif',
                                action="store_true",
                                help="produce animated gifs for graphs across the sampling factor range")
    parser_plot.add_argument('-l',
                                '--gene',
                                required=False,
                                type=str,
                                help="Provide gene-level csv for plotting")
    parser_plot.set_defaults(func=Plot.main)

    args=parser.parse_args()
    args.func(args)

def exit_handler(signum, frame):
    signal.signal(signal.SIGINT, original_sigint)

    try:
        while not len(childPIDS) == 0:
            for job in childPIDS:
                job.terminate()
                childPIDS.remove(job)
                print("quitting job")

    except KeyboardInterrupt:
        print("Ok ok, quitting")
        sys.exit(1)

    signal.signal(signal.SIGINT, exit_handler)

if __name__=="__main__":
    eva(sys.argv[1:])
