# PhacoBioIT: A scalable workflow solution for analysis and visualization ultra-low coverage NGS data from Pre-implantation genetic testing for aneuploidy.
## Introduction
PhacoBioIT was written by molecular biologist/bioinformatics expert anegin97. PhacoBioIT was created with the hope that other colleagues who are not specialized in using applications that require command line knowledge can perform the PGT-A data analysis process automatically with just the click of a button. The software is deployed as a web application, integrating data visualization tools along with a database to manage sample files.

PhacoBioIT mainly built with python and R, it includes the followings dependencies:
* Flask
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


**Clone the repository**
```bash
git clone https://github.com/Anegin24/PhacoBioIT.git
cd PhacoBioIT
```
**Prepare database and index:**
Download from link: ... 
```bash
unzip Ref
```
**Setup docker image:**
```bash
sudo docker build -t cnvviso -f CNVvisodocker .
sudo docker build -t phacobioit -f dockerfile .
```
**Run docker compose:**

In docker-compose.yml, you will see:

```bash "directory/DATA:/media/anegin97/DATA/DATA" ``` => Change "directory" to path contain DATA hard disk where store data
      
```bash "directory/PhacoBioIT/CNVdat:/home/anegin97/Bioinformatics/appforclient/CNVdat" ``` => Change "directory" to path contain DATA hard disk where store Rsqlite database
🚀🚀🚀🚀🚀
```bash
sudo docker compose up -d
```
