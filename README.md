## fasta2nuccomp

Summarize the contents of FASTA files.
This project consists of three main files.


- **fasta2nuccomp.py** a python script to summarize a FASTA file
- **fasta2nuccomp.Rmd** an RMarkdown script to create a report including visualization of the FASTA file summary
- **fasta2nuccomp.Rproj** an RStudio project to edit and knit the RMarkdown


The python script fasta2nuccomp.py takes a FASTA file as an argument and produces a per sequence tabular summary of the file.
The RMarkdown file fasta2nuccomp.Rmd demonstrates how to call the python script from R, calculate summary statistics, and create publication quality graphical summaries of the data.
Together these scripts demonstrate how to rapidly summarize FASTA files such as genome sequences.


---


![The yeast genome summarized with fasta2nuccomp.](S288C_genome.png)


