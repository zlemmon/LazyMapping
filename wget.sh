#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=1G

wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/lz_2_2_m_ACAGTG_L005_R1_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/lz_2_2_m_ACAGTG_L005_R2_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/Iz10m_1615_CTTGTA_L005_R1_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/Iz10m_1615_CTTGTA_L005_R2_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/lz_2_2_wt_TGACCA_L005_R1_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/lz_2_2_wt_TGACCA_L005_R2_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/Iz8m_1620_GCCAAT_L005_R1_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/Iz8m_1620_GCCAAT_L005_R2_001.fastq.gz
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/20161110_WKU_ParkSoonju_md5.txt
wget -q -c http://'user1':'!user!'@211.174.205.46/rawdata/20161110_WKU_ParkSoonju/20161110_WKU_ParkSoonju_Summary.txt
