# Count depth (bams are not shown in this example folder)
sh ../CollectDepthInfoForMSI.sh 20CF30397F.rmdup.bam 20CF30397F.info
sh ../CollectDepthInfoForMSI.sh 20CF30397B.rmdup.bam 20CF30397B.info
# Judge MSI
python ../JudgeMSI.py 20CF30397F.info 20CF30397B.info 20CF30397.MSI.csv