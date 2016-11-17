
# RNA-seq Analysis with Galaxy - Aligning Reads to a Bacteriophage Genome

The objective of this tutorial is to align a set of RNA-seq reads from the *Mycobacterium* phage Giles to its genome and perform a basic analysis (i.e., obtain a table of counts).

The overall process will be very similar to the [this guide](http://rpubs.com/kylescotshank/216067), which shows how to align the RNA-seq reads from the phage to its host's genome. The basic steps are as follows: 

  1. Import the reference genome and reference genome annotation.
  2. Import the RNA-seq reads.
  3. Perform diagnostic analyses of the RNA-seq reads.
  4. Align the RNA-seq reads to the reference genome.
  5. Generate a count of reads per gene that can be analyzed downstream with `R`.

__If you'd like to be able to see the output of the commands below in `R`, please follow the link above to the RPubs page.__

***

## Step 1: Import your Genome Information

In this step, we're going to download the genome assembly (a [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file) and the annotation (a [gff](http://useast.ensembl.org/info/website/upload/gff.html) file) onto our local machine from [EnsemblBacteria](http://bacteria.ensembl.org/index.html).

***

### Download the Genome Assembly

Click on [this link](http://applbio.mdibl.org/giles_phage.fa) to download the Giles genome assembly.

More information about the Giles genome assembly can be found [here](http://phagesdb.org/phages/Giles/)

***

### Download the Genome Annotation

Genome annotation can be one of the more difficult problems tackled in bioinformatics, mostly due to the plethora of file formats and transformation tools that are available. To simplify the task in this particular study, we have provided a suitable `gff` file, availble [here](http://applbio.mdibl.org/giles_phage.gtf).

***

### Import into Galaxy

From the Galaxy homepage, select __Get Data__.

![](Galaxy_GetData.tiff)

<br>
Next, click __Upload File from your computer__
  
   
![](Galaxy_UploadFileFromYourComputer.tiff)

<br>

Drag and drop both files that you wish to upload. Under `Type`, make sure to change the set the Geneome Assembly file to `fasta` and the Genome Annotation file to `gtf`. Click __Start__.

![](Galaxy_UploadScreen.tiff)

<br>
Click __close__. You can now see in your history bar (on the right) that you've successfully uploaded both files. Note that each item in your history equates to a "step" in your workflow - and each step is assigned a number in ascending order. This will be helpful downstream, as you will refer *back* to numeric steps in your workflow to feed data into the pipeline. 
***

## Step 2: Import RNA-seq Reads. 

The reads for the Giles RNA-Seq study were initially deposited in the [NCBI Short Read Archive](https://www.ncbi.nlm.nih.gov/sra).  The [European Nucleotide Archive (ENA)](http://www.ebi.ac.uk/ena) has mirrored those data and has made it easy to upload the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files into Galaxy.  

***

### Find the RNA-seq Reads.

Navigate to the ENA and enter the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) accession for this study, GSE43434. Click __Search__.

![](ENA_Search.tiff)

<br>
On the results page, select the first link (`SRP017906`) to look at the study. 

![](ENA_Results.tiff)

***

### Import into Galaxy

On the results page, make note of the `SRR-` values in the _Run accession_ column. Each of these files is an individual run through the sequencing machine. Copy the first string you see (`SRR647673`). 

![](ENA_ResultsList.tiff)

<br>
Return to Galaxy. From the toolbar, select __NCBI SRA Tools__. Then select __Extract reads in FASTQ/A format from NCBI SRA__.

![](Galaxy_NCBITools.tiff)

<br>
On this page, paste your copied run accession code (`SRR647673`) into the appropriate blank field. Then click execute. 

![](Galaxy_ExtractReads.tiff)

<br>
This will add the file (as a pending job) to your history on the right. Not that this process can take quite a bit of time to complete. Repeat the above for the other three `SRR-` files (`SRR647674`,`SRR647675`). Note that these files will be in the `fastqsanger` format.

***

## Step 3: Perform QC 

***

Quality control is an important step in bioinformatic (and all general data analysis) pipelines. We will be using `FastQC`. FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are:

  * Import of data from BAM, SAM or FastQ files (any variant)
  * Providing a quick overview to tell you in which areas there may be problems
  * Summary graphs and tables to quickly assess your data
  * Export of results to an HTML based permanent report
  * Offline operation to allow automated generation of reports without running the interactive application
  
***

### Running FastQC

From the Galaxy toolbar, click __NGS: QC and manipulation__. Then click __FastQC__.

![](Galaxy_FindFastQC.tiff)

<br>

Select one of your `FASTQ` files from the history to read in. Note that these files may have different referneces in your own history (in the example, the 3 `FASTQ` files imported from ENA are called 9: Extract Reads, 10: Extract Reads, and 11: Extra Reads). Then click __Execute__

![](Galaxy_RunFastQC.tiff)

***

### Examine FastQC Output

You will see two new additions to your history bar: a FastQC "RawData" file, and a FastQC "Webpage". You can download this file to your desktop and examine the FastQC Output. The main file of interest is the `html` file. Upon opening, you should see something similar to this:

![](FastQC.tiff)


***

## Step 4: Map to the Genome

Mapping refers to the process of aligning short reads to a reference sequence, whether the reference is a complete genome, transcriptome, or de novo assembly. There are numerous programs that have been developed to map reads to a reference sequence that vary in their algorithms and therefore speed. The program that we utilize in this pipeline is called `bowtie`. More information available [here](http://bowtie-bio.sourceforge.net/index.shtml)

***

### Running Bowtie

Click __NGS: Mapping__. Then click __Map with Bowtie for Illumina__.

![](Galaxy_FindBowtie.tiff)

<br>

In the first blank area ("_Will you select a reference genome..._?"), select __Use one from the history__. Then, choose your reference genome (`Giles_phage.fa`). Then, in the _FASTQ_ area, use one of your extracted reads from the previous steps. When you're ready, click __Execute__. 

![](Galaxy_UseBowtie.tiff)

<br>

Repeat this step for the remaining two `FASTQ` files. Note that your output files will be `SAM` files. For more information on `SAM`/`BAM` files, click [here](http://samtools.github.io/hts-specs/SAMv1.pdf)

***

## Step 5: Generate Counts per Read

To perform differential analysis, it's necessary to be able to calculate the number of reads mapping to each feature. Here, we think of a feature as an interval (i.e., a range of positions) on a chromosome or a union of such intervals. In the case of RNA-Seq, the features are typically genes, where each gene is considered here as the union of all its exons. One may also consider each exon as a feature, e.g., in order to check for alternative splicing. 

To perform this task, we will use the `htseq-count` program.

***

### Running htseq-count

From the toolbar, click __NGS: RNA Analysis__. Then click __htseq-count__.

![](Galaxy_FindHTSeq.tiff)

<br>

In the first input ("Aligned SAM/BAM file"), select one of the `SAM` files generated from `bowtie`. Then, make sure your `gtf` file is in the second input ("GFF File"). Leave the rest of the parameters as they are. Click __Execute__.

![](Galaxy_RunHTSeq.tiff)

When you've completed the above for all of your samples, don't forget to download the output files (which will be `.txt` files) to your local machine for downstream analysis!

***
