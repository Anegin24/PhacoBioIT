#!/bin/bash
cd /home/anegin97/Bioinformatics/appforclient/samplesheet
R -e "shiny::runApp('samplesheet.R',port = 2497,host='0.0.0.0')"
