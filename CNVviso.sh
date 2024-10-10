#!/bin/bash
cd /home/anegin97/Bioinformatics/appforclient/CNVviso
R -e "shiny::runApp('Vernewest240419.R',port = 2439,host='0.0.0.0')"
