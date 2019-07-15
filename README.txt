{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf400
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww21200\viewh11580\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #\
# photoferrotrophy.R\
1-D water column model for competition between oxygenic and anoxygenic phototrophs\
#\
This script models competition between oxygenic (cyanobacterial) and anoxygenic (photoferrotrophic) photosynthetic bacteria in an idealized 1-D water column.  The model parameters and equations are described in K. Ozaki, K.J. Thompson, R.L. Simister, S.A. Crowe, and C.T. Reinhard (2019) "Anoxygenic photosynthesis and the delayed oxygenation of Earth's atmosphere", Nature Communications, doi:10.1038/s41467-019-10872-z\
#\
As a courtesy, we request that others who use this code please cite Ozaki et al. (2019).  We also request that those who use/modify the code please send publications and/or modified code to the corresponding author (chris.reinhard@eas.gatech.edu)\
#\
REQUIREMENTS: R and/or R Studio, packages ReacTran, rootSolve, deSolve, shape, marelac, tictoc\
#\
TO RUN THE CODE\
From the command line:\
(1) navigate to the directory hosting photoferrotrophy.R\
(2) execute: Rscript photoferrotrophy.R\
(3) output will be stored in the as "photosynth.out.csv" in the directory specified in Line 200\
#\
Using R or R Studio:\
(1) open new console window\
(2) execute: source("~/\{DIRECTORY\}/photoferrotrophy.R")\
(3) output will be stored as above\
#}