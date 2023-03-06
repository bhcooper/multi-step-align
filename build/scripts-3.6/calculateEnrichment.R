#!/usr/bin/env Rscript

options(java.parameters = "-Xmx32G")
workDir = tempdir()

library(SELEX)
library(stringr)

selex.config(workingDir=workDir, maxThreadNumber=16)

R0.file = commandArgs(TRUE)[1]
data.file = commandArgs(TRUE)[2]
round = strtoi(commandArgs(TRUE)[3])
kLen = strtoi(commandArgs(TRUE)[4])

varLen = nchar(readLines(gzfile(data.file), 2)[2])

selex.defineSample('R0', R0.file, 'R0', 0, varLen, '', '')
selex.defineSample('RN', data.file, 'RN', round, varLen, '', '')

r0 = selex.sample(seqName = 'R0', sampleName = 'R0', round = 0)
r0.split = selex.split(r0)
r0.train = r0.split$train
r0.test = r0.split$test
dataSample = selex.sample(seqName = 'RN', sampleName = 'RN', round = round)

# MARKOV MODEL BUILT
kmax = selex.kmax(sample = r0.test, threshold = 100)
mm = selex.mm(sample = r0.train, order = NA, crossValidationSample = r0.test, Kmax = kmax, mmMethod="DIVISION")
mmscores = selex.mmSummary(sample = r0.train)
ido = which(mmscores$R==max(mmscores$R))
mm.order = mmscores$Order[ido]

selex.mmSummary()

# INFOGAIN USED TO CALCULATE KLEN
libLen = as.numeric(as.character(selex.getAttributes(dataSample)$VariableRegionLength))
selex.infogain(sample = dataSample, k = c(kLen), markovModel = mm)
aff = selex.affinities(sample=dataSample, k=kLen, markovModel=mm, numSort=FALSE)

aff = aff[,c(1,5,2,3,4,6)]
aff = aff[order(-aff$Affinity),]

write.table(aff, paste0(substr(data.file, 1, str_locate(data.file, "\\.")[1]-1), "_k", kLen, ".tsv"), row.names=FALSE, sep='\t', quote=FALSE)