# EVA - Expression Variation Analysis

# Sample Animated plot showing how the number of spots in the alignment affects the accuracy of expression estimation
![](https://github.com/alevar/EVA/blob/master/figures/analysis100:250/gene/png/boxID.gif)
# A general plot showing the median and 2nd/3rd quartile groups
![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/boxSFpa.png)
# A scatter matrix showing relationships between different statistical features
![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/scatterMatrixSF.png)
# Comparrison of rankings of transcripts by expression levels
![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/tauSF.png)

EVA is a tool designed to facilitate better understanding of the effects sequencing depth and coverage have on the estimation of gene-level and transcript-level expression(TPM).

The end goal of running EVA is to facilitate better understanding of variation and artifacts in transcript and gene-level expression estimates as a function of the number of reads sequenced and particular research interests (such as but not limited to: interest in transcripts with low or high expressions levels). We hope this tool will aid researchers in estimating the costs for their RNA-seq experiments.

As an example EVA was run on the full GTEX dataset in order to infer precision and cost tradeoffs of lower sequencing coverage for human transcriptomes.

Workflow:
1. "eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r auto"
	- this command parses through the alignment logs and builds information about each tissue and individual sample in the dataset and the total number of reads and percent of reads aligned.
	- the constructed dataset is saved locally as the output of the command
	- the auto range parameter binned the data and outputed the bin with the greatest mean number of aligned reads and saved this subset of the original dataset locally as well as a box plot for the number of aligned reads
	- after close inspection of both dataframes we came to the following conclusion:
		- samples suggested by the auto range parameter had very low standard deviation
		- the mean of 68M reads was too low as compared to many other samples in the original dataset which were distributed a little further from their means when binned
	- We ran the same command with custom ranges ("eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r minbound:maxbound") to establish the best bin for the task.
		- in our case we've got 12 samples in the range 90M-96M
2. "eva assemble -i <paths to the samples as outputed by eva select> -r 0.1:1.0:0.1 -o ./test -a <path to the annotation> -e <path to the reference> -w 10 -f 12 -t 2 -o <path to the output directory> -l"
	- eva assemble uses samtools, hisat and stringtie to assemble input reads and estimate expression levels. Stringtie is run with the "-e" flag which tells stringtie to "only estimate the abundance of given reference transcripts"
	- Explanation (more info in eva assemble --help):
		-i - perform for each of the given alignments
		-r - uses samtools hidden -S option to randomly downsample each alignment starting at 0.1 (10% retention of original reads) to 1.0 (all original reads retained - no downsampling) in 0.1 (10%) increments
		-o - output directory
		-w - randomness of each downsampling is based on a randomly generated seed. this option specifies how many times each downsampling with a different random seed should be done
		-f - split eva assemble into 12 forks each carrying out the pipeline for one assigned input sample
		-t - the number of cpu threads each fork is allowed to use. This gets passed to stringtie and samtools.
		-l - additionally estimates gene-level expression and stores separately.
	- For each assembly "eva assemble" will also extract information about transcripts and genes and will produce a single dataset for all assemblies (for each input: for each downsampling: for each repetition). The database contains the following information:
		- TISSUE CHROMOSOME REFERENCE_ID TRANSCRIPT_NAME START END SAMPLE FACTOR COVERAGE TPM
3. "eva analyze -l <path to the gene-level output of "eva analyze"> --de -o <path to the output dir>"
   "eva analyze -i <path to the transcript-level output of "eva analyze"> --de -o <path to the output dir>"
   - "eva analyze" outputs several datasets groupped by transcipts/genes and by the downsampling factor. The statistics included in these datasets can be plotted using "eva plot".

- Unique ID (tissueName:chromosome:referenceID:transcript/geneID:startLocation-endLocation) is assigned for each unique gene/transcript in the input dataset.

- Optionally if "-c" option is provided, genes/transcripts will be filtered by the specified TPM range. Thus it is possible to only look at a subset of the genes/transcripts.

- For each entry in the input dataset expression estimate is divided by the expression estimate of the control assembly of full set of aligned paired-end reads. This value is stored as %TPM

- False negatives are counted as the number of genes/transcripts that had control TPM>0 but TPM=0 at downsampling.
- False positives are counted as the number of genes/transcripts that had control TPM=0 but TPM>0 at downsampling.
- falseNegativesFull - number of reads which were present in the assembly of the original alignment but are no longer detected after downsampling in any of the repetitions of that particular downsampling (complete loss).
- the number of true positives and true negatives are computed likewise

- The following are computed based on raw expression values(TPM) as well as %TPM values
	- mean, coefficient of variation, standard deviation, median, normality of samples in each group
- For each sampling factor genes/transcripts are first grouped by uniqueID. Since we ran 10 random downsamplings of aligned reads to each sampling factor, each group contains 10 samples. Mean, coefficient of variation, standard deviation, median and normality were computed for each group. The resulting dataframe is saved on disk.

- All genes/transcripts in the dataset are grouped by sampling factor.
- Recall and precision are computed from the counts of true positives, true negatives, false positives, false negatives
- quartiles(25,50,75) and interquartile range are computed from %TPM values
- low and high whiskers as well as extremes were computed as follows: (include latex figures)
- standard deviation and coefficient of variation are computed when grouping genes/transcripts by sampling factor as well.

- For each of the downsamplings means of 10 random samples of each gene/transcipt are ranked. Kendall's tau correlation coefficient is then computed for each ranking against the ranking of the control assembly of full set of aligned paired-end reads.

- Reverse of a student t-test is performed on the dataset in order to determine the minumum expression value(TPM) which would be statistically signifficant under conservative two-tailed p-value of 0.01.

3. "eva plot -f <path to genes grouped by sampling factor> -i <path to genes grouped by unique id> -g -o <path to the output dir>"
   "eva plot -f <path to transcripts grouped by sampling factor> -i <path to transcripts grouped by unique id> -g -o <path to the output dir>"

	- The following plots are produced by "eva analyze" for the data grouped by sampling factor:
			- tauSF - kendal's tau correlation coeeficient plots with a shared y-axis (tau coeeficient). Plots include:
				- all transcripts
				- top 10% of transcripts by percentAway value
				- top 20% of transcripts by percentAway value
				- top 50% of transcripts by percentAway value
				- bottom 50% of transcripts by percentAway value
				- bottom 20% of transcripts by percentAway value
				- bottom 10% of transcripts by percentAway value
			- boxSF - median and quartiles
			- PCA - a biplot of principal component 1 vs principal component 2. The biplot also shows the direction and weights of the vectors for each of the following input features:
				- falsePositives
				- falseNegatives
				- falseNegativesFull
				- median
				- weightedNormalizedNumExtremes
				- std
				- tauFull

        - the following plots are produced for each sampling factor:
        	- boxID - median and quartiles. X-axis: coverage; y-axis: percentAway
        	- groupedID - each statistical measure plotted separately against the coverage (as reported for unsampled assembly)
        	- scatterMatrixByID - scatter matrix of the measures against each other. Diagonal shows histograms
        	- Additionally the each measure was transformed via BoxCox transformation and the same plots were produced for each downsampling from the transformed data.

4. By analysing the results we came to the following conclusions:
	- To be updated in a few days. The final run is still working.

Please ensure the packages from the requirements.txt are installed in your python environment.
