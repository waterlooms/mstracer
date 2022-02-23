#!/bin/bash
"""
This script is to highlight the flow of running WaterlooMS from RAW data acquisition to mstracer output.
"""

mkdir data/PXD005573/

wget -P data/PXD005573/ ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005573/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.raw
wget -P data/PXD005573/ ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005573/uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta


# Convert the data set, with centroiding
docker run -it -v $PWD/data/PXD005573:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:x64 wine msconvert --zlib --filter "peakPicking true 1-" --mzML /data/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.raw

# Convert the file to be owned by the USER:GROUP of the host 
docker run -it -v $PWD/data/PXD005573:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:x64 chown $(id -u ${USER}):$(id -g ${USER}) Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.mzML


# Copy the filename for the data over for a more weildy filename (optional)
mv data/PXD005573/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.mzML data/PXD005573/dia_data.mzML

# Invoke WaterlooMS on the now converted mzMmL file
# Run the Maven lifecycle to generate a FatJar of WaterlooMSs
mvn clean install package



docker run --rm -v $PWD:/app -w /app trackerrr/mstracer java -jar /mstracer/mstracer.jar -mzXML [PATH_TO_MZXML]



"""
For DEV you want to mount the entire directory and have access to the file(s)
However for PROD you want only the /src/main/python directory mounted





"""








