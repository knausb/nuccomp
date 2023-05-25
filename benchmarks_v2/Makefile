
# https://kbroman.org/minimal_make/

# https://stackoverflow.com/a/49908054

knit :
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_Arabidopsis.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_cottonwood.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_hazelnut.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_hemp.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_mint.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_potato.Rmd')"
	Rscript -e "rmarkdown::render('fasta2nuccomp_v2_rice.Rmd')"


        
#FASTA_SOTU = GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#FASTA_ARTH = GCF_000001735.4_TAIR10.1_genomic.fna.gz


get_fa :
# Arabidopsis
ifneq ("$(wildcard GCF_000001735.4_TAIR10.1_genomic.fna.gz)","")
#	FILE_EXISTS = 1
	@echo "File 'GCF_000001735.4_TAIR10.1_genomic.fna.gz' already exists"
else
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz --output GCF_000001735.4_TAIR10.1_genomic.fna.gz
endif
# Cottonwood
ifneq ("$(wildcard GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz)","")
	@echo "File 'GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz' already exists"
else
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/775/GCF_000002775.5_P.trichocarpa_v4.1/GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/775/GCF_000002775.5_P.trichocarpa_v4.1/GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz --output GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz
endif
# Hazelnut
ifneq ("$(wildcard GCA_901000735.2_CavTom2PMs-1.0_genomic.fna.gz)","")
	@echo "File 'GCA_901000735.2_CavTom2PMs-1.0_genomic.fna.gz' already exists"
else
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/901/000/735/GCA_901000735.2_CavTom2PMs-1.0/GCA_901000735.2_CavTom2PMs-1.0_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/901/000/735/GCA_901000735.2_CavTom2PMs-1.0/GCA_901000735.2_CavTom2PMs-1.0_genomic.fna.gz --output GCA_901000735.2_CavTom2PMs-1.0_genomic.fna.gz
endif
# Hemp
ifneq ("$(wildcard GCF_900626175.2_cs10_genomic.fna.gz)","")
	@echo "File 'GCF_900626175.2_cs10_genomic.fna.gz' already exists"
else
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/626/175/GCF_900626175.2_cs10/GCF_900626175.2_cs10_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/626/175/GCF_900626175.2_cs10/GCF_900626175.2_cs10_genomic.fna.gz --output GCF_900626175.2_cs10_genomic.fna.gz
endif
# Mint
ifneq ("$(wildcard GCA_001642375.2_Mlong_CMEN585_v3_genomic.fna.gz)","")
	@echo "File 'GCA_001642375.2_Mlong_CMEN585_v3_genomic.fna.gz' already exists."
else
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/642/375/GCA_001642375.2_Mlong_CMEN585_v3/GCA_001642375.2_Mlong_CMEN585_v3_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/642/375/GCA_001642375.2_Mlong_CMEN585_v3/GCA_001642375.2_Mlong_CMEN585_v3_genomic.fna.gz --output GCA_001642375.2_Mlong_CMEN585_v3_genomic.fna.gz
endif
# Potato
ifneq ("$(wildcard GCF_000226075.1_SolTub_3.0_genomic.fna.gz)","")
	@echo "File 'GCF_000226075.1_SolTub_3.0_genomic.fna.gz' already exists"
else
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/226/075/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/226/075/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.fna.gz --output GCF_000226075.1_SolTub_3.0_genomic.fna.gz
endif
# Rice
ifneq ("$(wildcard GCF_001433935.1_IRGSP-1.0_genomic.fna.gz)","")
	@echo "File 'GCF_001433935.1_IRGSP-1.0_genomic.fna.gz' already exists"
else
	#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz --output GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
endif




clean :
#	@rm -f *.fna.gz
	@rm -f *_nuccomp.csv
	@rm -f *_CG_wins.bed
	@rm -f *_nuccomp_summary.csv
#	@rm -f *.html
#ifneq ("$(wildcard GCF_000226075.1_SolTub_3.0_genomic.fna.gz)","")
#	rm GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#	FILE_EXISTS = 1
#else
#	@echo "File 'GCF_000226075.1_SolTub_3.0_genomic.fna.gz' does not exist."
#	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/226/075/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.fna.gz
#	FILE_EXISTS = 0
#endif
#	rm GCA_019925245.1_MQ-UA-ANU_Oaus-KR_1.0_genomic.fna.gz
#	rm GCF_000001735.4_TAIR10.1_genomic.fna.gz



# EOF.
