# PhacoBioIT: A scalable workflow solution for analysis and visualization ultra-low coverage NGS data from Pre-implantation genetic testing for aneuploidy.
## Introduction
PhacoBioIT was written by molecular biologist/bioinformatics expert anegin97. PhacoBioIT was created with the hope that other colleagues who are not specialized in using applications that require command line knowledge can perform the PGT-A data analysis process automatically with just the click of a button. The software is deployed as a web application, integrating data visualization tools along with a database to manage sample files.

PhacoBioIT mainly built with python and R, it includes the followings dependencies:
* bwa
* samtools
* zip
* unzip
* cnvkit
* ggplot2
* tidyverse
* dplyr
* shiny
* DT
* RSQLite
* ggtext
* dplyr
* plotly
* smoother
* magick
* patchwork
* purrr
* here
## Installation
### Docker
The pipeline can be installed via Docker as well.
```bash
git clone https://github.com/Anegin24/PhacoBioIT.git
cd PhacoBioIT
sudo docker build -t PhacoBioIT .
```

The data input directory from the container is /media/anegin97/DATA/DATA, user can run PhacoBioIT commands by mounting the host database (where the git cloned into, e.g. ~/workspace/amromics) into this destination (by using -v).
```bash
chmod 777 /directory/DATA
sudo docker run -d -p 9999:9999 -p 80:80 --name PhacoBioIT  -v /media/anegin97/DATA/DATA:/media/anegin97/DATA/DATA PhacoBioIT
```
## Usage
