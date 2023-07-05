# Table of Contents
1. [nuccomp](#nuccomp)
2. [Installation](#installation)
3. [Input files](#input-files): FAST[AQ] format
4. [Use_case_1](#use-case-1): command line python
5. [Use case 2](#use-case-2): Calling python from Rmarkdown


## nuccomp

Summarize the contents of FAST[AQ] files.
This project consists of three main files.


- **nuccomp.py** a python script to summarize a FASTA file
- **nuccomp.Rmd** an RMarkdown script to create a report including visualization of the FASTA file summary
- **nuccomp.Rproj** an RStudio project to edit and knit the RMarkdown


The python script nuccomp.py takes a FASTA or FASTQ file as an argument and produces a per sequence tabular summary of the file.
The RMarkdown file nuccomp.Rmd demonstrates how to call the python script from R, calculate summary statistics, and create publication quality graphical summaries of the data.
Together these scripts demonstrate how to rapidly summarize FAST[AQ] files such as genome sequences.


---


![The yeast genome summarized with nuccomp.](S288C_reference_sequence_R64-2-1_20150113.png)



## Installation


In order to use this code the user will need python and R installed on their system.
We also recommend using the integrated development environment [RStudio](https://posit.co/products/open-source/rstudio/).
Please consult the respective project pages to ensure these are installed.

The python script uses [Biopython](https://biopython.org/).
Please ensure that this is installed.
The following is an example of what to expect if Biopython is not installed.


```
$ ./nuccomp.py 
Traceback (most recent call last):
  File "./nuccomp.py", line 6, in <module>
    from Bio import SeqIO
ModuleNotFoundError: No module named 'Bio'
```


Please see Biopython's documentation to make sure it is installed correctly.
Note that a system may include several versions of python (including conda environments).
This can create a situation where Biopython, or other dependencies, are installed in one environment but not in the environment a user is attempting to use at any particular time.
This means the user needs to ensure that dependencies are installed into the environment they are attempting to use.
If you are using conda you may want to include something like the following in your `~/.Rprofile` in order to help R find your desired version of python.


```
Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/biopython/bin/python")
```


Installation of `nuccomp` can be accomplished by cloning the GitHub repository as follows.


```
git clone git@github.com:knausb/fasta2nuccomp.git
```

Or downloading a release.
Once on a local machine the scripts should be placed in the user's path.


Each python script should return a usage message when called with no arguments, and a brief help message when the '-h' flag is included.

```
$ nuccomp.py
usage: nuccomp.py [-h] [-v] INFILE
nuccomp.py: error: the following arguments are required: INFILE
$ nuccomp.py -h
usage: nuccomp.py [-h] [-v] INFILE

Determine nucleotide composition of FASTA files. Requires python >=3.7.3 and
Biopython >= 1.78.

positional arguments:
  INFILE         FASTA file containing nucleotides.

options:
  -h, --help     show this help message and exit
  -v, --verbose  increase output verbosity
```



## Input files

Nucleotide data is typically stored in a FASTA format file.
This is a text file with a specific format.
A description line begins with '>' and is followed by one or more lines of sequence.
You can learn more about the FASTA format at Wikipedia's [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) page.
The yeast (*Saccharomyces cerevisiae*) strain S288C genome is included with ```nuccomp``` as an example FASTA file.
Below is a brief example of this file.
Note that the sequence portion of each record may exist on a single line or may be spread over many lines.
Because the sequence portion may be spread over many lines it is recommended that applications that are specifically designed to handle FAST files be used as opposed to treating these files as simple text.


```
>ref|NC_001133| [org=Saccharomyces cerevisiae]
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA
CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT
>ref|NC_001134| [org=Saccharomyces cerevisiae]
AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGATGTTCAACCA
AAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATAGGTTGTCTTTTTATCCCACT
```

While nuccomp.py can also be used to process FASTQ files.
For example, we use this script to summarize FASTQ sequencing results received from a sequencing center and prior to assembly.
This provides the opportunity to determine the read count and sequence length distribution in these sequencing libraries prior to assembly efforts.
The execution time for sequence_comp.py scale with file size, so processing sequencing libraries should take longer than processing assemblies, but our experience indicates this to be feasible.
Below is an example of a FASTQ file, more information can be found on Wikipedia's [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) page.
In the FASTQ format there are two description lines, prefixed with '@' on the first line and then '+' on the third line.
The sequence is on the second line, and the quality string is on the fourth line.


```
@SRR10238608.1 1 length=22221
AAATTGATGATTGGCCATTTTGTTTTCTAGAAGGTGCAACTGATCAGGATGAGGCTAGAGCTATTGTCTGATTGCTATAG
+SRR10238608.1 1 length=22221
Y{|e~mRfmVF_K~]~nF}|tjI}mgjutjIr:panhEkdpntsZmLp`Yj\[~idWgT?HbU>[JjYcnnThijfisvj
@SRR10238608.2 2 length=20840
TTATAAGGGTATTTATGTAAATCCCCTTTATATATTAGTCATAGGGTGTATGAGTTCGTATTAAAGACTATAAATAGACC
+SRR10238608.2 2 length=20840
~~~~~~~~~~~_~~~~~~R~~~|~~~R~~~z~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~e~~~~~~~~X~~~~~~~~
```


## Use case 1
**Command line python**

The python script `nuccomp.py` can be invoked, using a FASTA or FASTQ file that may be gzipped, as follows.


```
$ ./python/nuccomp.py S288C_reference_sequence_R64-2-1_20150113.fsa.gz
```

This example was intended to be executed from the root directory for `nuccomp`.
This should result in the following file in the current working directory.


```
S288C_reference_sequence_R64-2-1_20150113_nuccomp.csv
```

This is the input file name with the FAST[AQ] and any 'gz' extension removed and the suffix `_nuccomp.csv` added.
This file summarizes the contents of the FAST[AQ] file.


Processing of this file can be accomplished by following the example code in `nuccomp.html`.
Note that GitHub does not 'serve' these *.html files so the user will have to view a local copy.
Alternatively, the file `nuccomp.Rmd` can be compiled, or knit, into *.html using RStudio, or from the command line as follows. 

```
Rscript -e "rmarkdown::render('nuccomp.Rmd')"
```


## Use case 2
**Calling python from Rmarkdown**

R can be integrated with python allowing the python script to be run as a part of the Rmarkdown compilation.
This may require configuration that is beyond the scope of this project.
It has been our experience that installing the `reticulate` package will be required.

```
install.packages('reticulate')
```

Adding the below line to your `~/.Rprofile` file can help R know which version of python to use.

```
Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/biopython/bin/python")
```


Once R and python are configured the `nuccomp.Rmd` file can be modified to evaluate the code chucks containing python code.
This will allow the python code to be called from the Rmarkdown allowing processing to occur in one step.


