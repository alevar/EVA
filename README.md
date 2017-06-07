# EVA - Expression Variation Analysis

EVA is a tool designed to facilitate better understanding of the effects sequencing depth and coverage have on the estimation of gene-level and transcript-level expression(TPM).

The end goal of running EVA is to facilitate better understanding of variation and artifacts in transcript and gene-level expression estimates as a function of the number of reads sequenced and particular research interests (such as but not limited to: interest in transcripts with low or high expressions levels). We hope this tool will aid researchers in estimating the costs for their RNA-seq experiments.

As an example EVA was run on the full GTEX dataset in order to infer precision and cost tradeoffs of lower sequencing coverage for human transcriptomes.

![](https://github.com/alevar/EVA/blob/master/figures/pubFOLD_RECALL.png)
###Figure. Plot shows the change in variation of estimated gene expression levels (%TPM) as a function of the number of aligned reads used in the assembly. Variation in gene expression estimates is measured as percent from the control assembly of full set of aligned paired-end reads. The plot shows the median (dark line), second and third quartiles (dark-red area) and whiskers (light-red area).

![](https://github.com/alevar/EVA/blob/master/figures/pubVARIATION.png)
###Figure. Plot shows the change in recall (blue) and number of genes with significant change in expression estimation (green) as a function of the number of aligned reads used in the assembly. Genes counted for fold change graph were selected as those having two or more fold increase or decrease in TPM from the control assembly of full set of aligned paired-end reads.

##Please ensure the packages from the requirements.txt are installed in your python environment.

#Workflow:
1 "eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r auto"
	- this command parses through the alignment logs and builds information about each tissue and individual sample in the dataset and the total number of reads and percent of reads aligned.
	- the constructed dataset is saved locally as the output of the command
	- the auto range parameter binned the data and outputed the bin with the greatest mean number of aligned reads and saved this subset of the original dataset locally as well as a box plot for the number of aligned reads
	- after close inspection of both dataframes we came to the following conclusion:
		- samples suggested by the auto range parameter had very low standard deviation
		- the mean of 68M reads was too low as compared to many other samples in the original dataset which were distributed a little further from their means when binned
	- We ran the same command with custom ranges ("eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r minbound:maxbound") to establish the best bin for the task.
		- in our case we've got 12 samples in the range 90M-96M
2 "eva assemble -i <paths to the samples as outputed by eva select> -r 0.1:1.0:0.1 -o ./test -a <path to the annotation> -e <path to the reference> -w 10 -f 12 -t 2 -o <path to the output directory> -l"
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
3 "eva analyze -l <path to the gene-level output of "eva analyze"> --de -o <path to the output dir>"
   "eva analyze -i <path to the transcript-level output of "eva analyze"> --de -o <path to the output dir>"
   - "eva analyze" outputs several datasets grouped by transcipts/genes and by the downsampling factor. The statistics included in these datasets can be plotted using "eva plot".

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

4 "eva plot -f <path to genes grouped by sampling factor> -i <path to genes grouped by unique id> -g -o <path to the output dir>"
   "eva plot -f <path to transcripts grouped by sampling factor> -i <path to transcripts grouped by unique id> -g -o <path to the output dir>"
	- The following plots are produced by "eva analyze" for the data grouped by sampling factor:
			## Kendal's tau ranking correlations
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/tauSF.png)
			###Figure. Plot shows the change in gene ranking by expression estimate (TPM) disagreement from ranking of the control assembly produced from a full set of aligned paired-end reads.
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/transcript/png/tauSF.png)
			###Figure. Plot shows the change in transcript ranking by expression estimate (TPM) disagreement from ranking of the control assembly produced from a full set of aligned paired-end reads.

			## Change in expression variation(%TPM)
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/boxSFpa.png)
			###Figure. Plot shows the change in variation of estimated gene expression levels (%TPM) as a function of the number of aligned reads used in the assembly. Variation in gene expression estimates is measured as percent from the control assembly of full set of aligned paired-end reads. The plot shows the median (dark line), second and third quartiles (dark-red area) and whiskers (light-red area).
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/transcript/png/boxSFpa.png)
			###Figure. Plot shows the change in variation of estimated transcript expression levels (%TPM) as a function of the number of aligned reads used in the assembly. Variation in transcript expression estimates is measured as percent from the control assembly of full set of aligned paired-end reads. The plot shows the median (dark line), second and third quartiles (dark-red area) and whiskers (light-red area).

			## Significant change in expression estimates
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/foldIncrease.png)
			###Figure. Plot shows the number of genes with significant change in expression estimation as a function of the number of aligned reads used in the assembly. Genes were selected as those having two or more fold increase or decrease in TPM from the control assembly of full set of aligned paired-end reads.
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/transcript/png/foldIncrease.png)
			###Figure. Plot shows the number of transcripts with significant change in expression estimation as a function of the number of aligned reads used in the assembly. Genes were selected as those having two or more fold increase or decrease in TPM from the control assembly of full set of aligned paired-end reads.

			## Recall vs Precision
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/gene/png/recallPrecision.png)
			###Figure. Plot shows the change in fraction of expressed genes recovered in assembly (Recall) as a function of the number of aligned reads used in the assembly.
			![](https://github.com/alevar/EVA/blob/master/figures/analysisFull/transcript/png/recallPrecision.png)
			###Figure. Plot shows the change in fraction of expressed transcripts recovered in assembly (Recall) as a function of the number of aligned reads used in the assembly.

			Some other plots produced by "eva plot" are:
			- pca1 vs pca2 biplot
			- normality of random samples in each unique ID group
			- scattermatrix of all statistics computed from grouped data
			- falseNegatives
			- Contiguous boxplot based on coefficient of variation

        - Additionally several plots are produced for each sampling factor based on the data grouped by unique gene/transcript ID:
        	## Change in individual gene expression estimation
			![](https://github.com/alevar/EVA/blob/master/figures/analysis100:250/gene/png/boxID.gif)
			###Figure. Animated plot shows the change in individual gene expression estimate variation as a function of the number of paired-end reads in the alignment. Only genes with control expression estimate value (TPM) greater than 100 and smaller than 250 were used to produce this plot.
			![](https://github.com/alevar/EVA/blob/master/figures/analysis100:250/transcript/png/boxID.gif)
			###Figure. Animated plot shows the change in individual transcript expression estimate variation as a function of the number of paired-end reads in the alignment. Only transcripts with control expression estimate value (TPM) greater than 100 and smaller than 250 were used to produce this plot.

        	- groupedID - each statistical measure plotted separately against the coverage (as reported for unsampled assembly)
        	- scatterMatrixByID - scatter matrix of the measures against each other. Diagonal shows histograms
        	- Additionally the each measure was transformed via BoxCox transformation and the same plots were produced for each downsampling from the transformed data.