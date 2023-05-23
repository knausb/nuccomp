
# https://kbroman.org/minimal_make/

# https://stackoverflow.com/a/49908054

knit :
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_Arabidopsis.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_oryza.Rmd')"

objects = GCF_000226075.1_SolTub_3.0_genomic.fna.gz GCF_000001735.4_TAIR10.1_genomic.fna.gz

list : $(objects)/*
	for file in $^ ; do \
		echo "Hello" $${file} ; \
	done
        
#FASTA_SOTU = GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#FASTA_ARTH = GCF_000001735.4_TAIR10.1_genomic.fna.gz


get_fa :
#	objects = GCF_000226075.1_SolTub_3.0_genomic.fna.gz
ifneq ("$(wildcard GCF_000226075.1_SolTub_3.0_genomic.fna.gz)","")
#	FILE_EXISTS = 1
	@echo "File 'GCF_000226075.1_SolTub_3.0_genomic.fna.gz' already exists"
	@echo "File 'GCF_000226075.1_SolTub_3.0_genomic.fna.gz' already exists $(objects)"
else
	@echo "File does not exist. objects"
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/226/075/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#	FILE_EXISTS = 0
endif


#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/925/245/GCA_019925245.1_MQ-UA-ANU_Oaus-KR_1.0/GCA_019925245.1_MQ-UA-ANU_Oaus-KR_1.0_genomic.fna.gz
#ifneq ("GCF_000001735.4_TAIR10.1_genomic.fna.gz")
#ifneq ( "$(wildcard GCF_000001735.4_TAIR10.1_genomic.fna.gz)", )
#	FILE_EXISTS = 1
#else
#	FILE_EXISTS = 0
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
#endif






clean :
ifneq ("$(wildcard GCF_000226075.1_SolTub_3.0_genomic.fna.gz)","")
	rm GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#	FILE_EXISTS = 1
else
	@echo "File 'GCF_000226075.1_SolTub_3.0_genomic.fna.gz' does not exist."
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/226/075/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#	FILE_EXISTS = 0
endif
#	rm GCA_019925245.1_MQ-UA-ANU_Oaus-KR_1.0_genomic.fna.gz
#	rm GCF_000001735.4_TAIR10.1_genomic.fna.gz



# EOF.
