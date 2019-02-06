#######################
#gating_Double_Fluorescent_Strains.R       
#
#started: 01/07/2016
#modified: 2/6/2019
#
#author1: G Avecilla, S Lauer, D Gresham
#author2: N Brandt
######################

######################
#This script is a generic file specific for analyzing the data from flow cytometry in the FCS3 formate
#For your exact purposes, the script must be modified accordingly. 

###########################################################################################################################
#This script is intended to read in .fcs files and perform manual gating for i) single cells, ii) debris and iii) fluorescence 
#The fcs files should be in a folder with a name descriptive of the experiment and have a assocated sample sheet with the name format of samplesheet_<name>
#Should include at least the following:
#       * column1 = Well
#       * column2 = Strain

#
#Gating is performed with untransformed data
#
#Individual gates are saved in a file gates.Rdata for use with the Gresham Lab Flow Cytometry Analysis.Rmd pipeline
###########################################################################################################################

##To be run the first time if packages are not installed.
#source("http://bioconductor.org/biocLite.R")
#biocLite("flowViz")
#biocLite("flowCore")
#biocLite("ggcyto")

#Load libraries
library(flowCore)
library(flowViz)
library(ggcyto)
library(ggforce)

#Read in the data
#Set working directory to the folder in which you have stored your .fcs files
#Read in all the fcs files in the directory.

#working directory
dir = '.'

#file location
path.data = paste(getwd(),"/",sep="")

#set name of run to create gates for
name <- "Aria_Run"

#load sample sheet
sample.sheet <- read.csv(paste(path.data,"samplesheet_",name,".csv", sep=""))

#read in fcs files in order presented in sample sheet (based on well identifier)
files <- paste(path.data,name,"/",sort(factor(list.files(paste(path.data,name,"/", sep=""),full.names=FALSE), levels = paste(sample.sheet$Well,".fcs",sep="" ), ordered=TRUE)),sep="")
flowData <- read.ncdfFlowSet(files=files, pattern=".fcs", alter.names = TRUE)


#rename sample name of flow set to make it easier to identify
sampleNames(flowData) <- paste(gsub(" ","_",sample.sheet$Strain),"_",sub(" ","_",sample.sheet$Well), sep="")

#set controls on which to gate on
#these are some of the base types of controls you may want for each channel you are taking measurements from 
neg.signal <- 1
mCitrineGAP1.signal <- 2
mCherryMEP2.signal <- 3
mCherryHXT6_7.signal <- 6

##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area
#******Please note you may need to adjust the x and y plot limits to properly visualize your data

plot(flowData[[neg.signal]], c('FSC.H','FSC.A'), xlim=c(0,3e6), ylim=c(0,3e6),smooth=T)

singlet.gate <- locator(100, type='l', col='red')
gm.1 <- matrix(,length(singlet.gate$x),2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- singlet.gate$x
gm.1[,2] <- singlet.gate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

ggcyto(flowData[neg.signal], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + geom_gate(pg.singlets)

#Look at the gating on the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#test that the singlet gate looks reasonable for All samples
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)


##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.
#******Please note you may need to adjust the x and y plot limits to properly visualize your data


plot(flowData[[neg.signal]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,1e6),smooth=T)

debris.gate <- locator(100, type='l', col='red')
gm.2 <- matrix(,length(debris.gate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- debris.gate$x
gm.2[,2] <- debris.gate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#Look at the gating on the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)


#test that the singlet gate looks reasonable for All samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)



#############################
#3. FLUORESCENCE####
#******Please note you may need to adjust the x and y plot limits to properly visualize your data
# You may need to add or modify these gates depending on how you are analyzing the signal and how many channels you are working with
#This is the base code block you can use to add more gates CONTROLVARIABLE, TYPEOFSIGNAL, GATENUMBERINCODE, TYPEOFSIGNAL, should all be replaced with vairables that make sense in your code
######START
##Plot 
#     plot(flowData[[CONTROLVARIABLE]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e4),smooth=T)
##draw the gate and create a gate ploygon
#      TYPEOFSIGNAL.gate <- locator(100, type='l', col='red')
#      gm.GATENUMBERINCODE <- matrix(,length(TYPEOFSIGNAL.gate$x),2)
#      colnames(gm.GATENUMBERINCODE) <- c('FSC.A','FL1.A')
#      gm.GATENUMBERINCODE[,1] <- TYPEOFSIGNAL.gate$x
#      gm.GATENUMBERINCODE[,2] <- TYPEOFSIGNAL.gate$y
#      gate.TYPEOFSIGNAL <- polygonGate(filterId="TYPEOFSIGNAL",.gate=gm.GATENUMBERINCODE)
##Overlay and check the new gate
#     ggcyto(flowData[c(CONTROLVARIABLE)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.TYPEOFSIGNAL)
######END

##Plot the control sample that has non-fluorescing cells in FL1 channel
plot(flowData[[negFL1.signal]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e4),smooth=T)

negFL1signal.gate <- locator(100, type='l', col='red')
gm.3 <- matrix(,length(negFL1signal.gate$x),2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- negFL1signal.gate$x
gm.3[,2] <- negFL1signal.gate$y
gate.negFL1 <- polygonGate(filterId="zerosignal",.gate=gm.3)

##Overlay and check the new gate
ggcyto(flowData[c(neg.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1)

#Look at the gating on all the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)



##Draw a new gate for gating gate for mCitrineGAP1 flourescent cells
plot(flowData[[mCitrineGAP1.signal]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e5),smooth=T)
polygon(negFL1signal.gate)

mCitrineGAP1.gate <- locator(100, type='l', col='blue')
gm.4 <- matrix(,length(mCitrineGAP1.gate$x),2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,1] <- mCitrineGAP1.gate$x
gm.4[,2] <- mCitrineGAP1.gate$y
gate.mCitrineGAP1 <- polygonGate(filterId="mCitrineGAP1_signal",.gate=gm.4)

##Overlay and check the new gate with the old gate
ggcyto(flowData[mCitrineGAP1.signal], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1)

#Look at the gating on all the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.negFL1) + geom_gate(gate.mCitrineGAP1) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)


##Plot the control sample that has non-fluorescing cells in FL2 channel
plot(flowData[[neg.signal]], c('FSC.A','FL2.A'), xlim=c(0,3e6), ylim=c(0,5e3),smooth=T)

negFL2signal.gate <- locator(100, type='l', col='red')
gm.3 <- matrix(,length(negFL2signal.gate$x),2)
colnames(gm.3) <- c('FSC.A','FL2.A')
gm.3[,1] <- negFL2signal.gate$x
gm.3[,2] <- negFL2signal.gate$y
gate.negFL2 <- polygonGate(filterId="zerosignal",.gate=gm.3)

##Overlay and check the new gate
ggcyto(flowData[c(neg.signal)], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2)

#Look at the gating on all the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e3) + geom_gate(gate.negFL2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)


##Draw a new gate for gating gate for mCherryMEP2 flourescent cells
plot(flowData[[mCherryMEP2.signal]], c('FSC.A','FL2.A'), xlim=c(0,3e6), ylim=c(0,1e4),smooth=T)
polygon(negFL2signal.gate)

mCherryMEP2.gate <- locator(100, type='l', col='blue')
gm.4 <- matrix(,length(mCherryMEP2.gate$x),2)
colnames(gm.4) <- c('FSC.A','FL2.A')
gm.4[,1] <- mCherryMEP2.gate$x
gm.4[,2] <- mCherryMEP2.gate$y
gate.mCherryMEP2 <- polygonGate(filterId="mCherryMEP2_signal",.gate=gm.4)

##Overlay and check the new gate with the old gates
ggcyto(flowData[mCherryMEP2.signal], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) 

#Look at the gating on all the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryMEP2) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)

##Draw a new gate for gating gate for mCherryHXT6/7 flourescent cells
plot(flowData[[mCherryHXT6_7.signal]], c('FSC.A','FL2.A'), xlim=c(0,3e6), ylim=c(0,1e4),smooth=T)
polygon(negFL2signal.gate)

mCherryHXT6_7.gate <- locator(100, type='l', col='blue')
gm.4 <- matrix(,length(mCherryHXT6_7.gate$x),2)
colnames(gm.4) <- c('FSC.A','FL2.A')
gm.4[,1] <- mCherryHXT6_7.gate$x
gm.4[,2] <- mCherryHXT6_7.gate$y
gate.mCherryHXT6_7 <- polygonGate(filterId="mCherryHXT6_7_signal",.gate=gm.4)

##Overlay and check the new gate with the old gates
ggcyto(flowData[mCherryMEP2.signal], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) 

#Look at the gating on all the controls
ggcyto(flowData[c(neg.signal,mCitrineGAP1.signal,mCherryMEP2.signal,mCherryHXT6_7.signal)], aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 2)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 3)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 4)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 5)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 6)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 7)
ggcyto(flowData, aes(x = `FSC.A`, y =  `FL2.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e4) + geom_gate(gate.negFL2) + geom_gate(gate.mCherryHXT6_7) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 8)



#Save the gate information to an R data file
#Remeber to add any additional gate variables to your save file
rm(list=c("flowData")) 
save(singlet.gate, pg.singlets, debris.gate, pg.nondebris, negFL1signal.gate, gate.negFL1, mCitrineGAP1.gate, gate.mCitrineGAP1, negFL2signal.gate, gate.negFL2, mCherryMEP2.gate, gate.mCherryMEP2, mCherryHXT6_7.gate, gate.mCherryHXT6_7, file=paste(name,"_gates_",Sys.Date(),".Rdata",sep=""))
