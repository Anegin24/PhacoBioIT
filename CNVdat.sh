#!/bin/bash
cd /home/anegin97/Bioinformatics/appforclient/CNVdat
R -e "shiny::runApp('CNVdat.R',port = 2437,host='0.0.0.0')"
