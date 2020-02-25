mkdir -p results
seqsimla -popfile EUR_chr1.bed.gz -recfile EUR_chr1.rec -famfile simped1000.ped -folder results -header Sim1 -batch 100 -site 80 --mode-prev -prev 0.1 -or 2.0
