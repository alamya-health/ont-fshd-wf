mkdir -p ./OC195M_FSHD_debug
cd ./OC195M_FSHD_debug

mkdir -p /tmp/OC195M_FSHD_debug
cd /tmp/OC195M_FSHD_debug

aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/build-fshd-report/OC195M/OC195M.fshd.report.html" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/build-fshd-report/OC195M/OC195M.fshd.report.summary.tsv" .

aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/classify-fshd-reads/OC195M/OC195M.fshd.classification/A_4qA_all-reads.csv" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/classify-fshd-reads/OC195M/OC195M.fshd.classification/FSHD_overview-statistics.csv" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/classify-fshd-reads/OC195M/OC195M.fshd.classification/blast_results/4qA_reads_blast.csv" .

aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/profile-fshd-methylation/OC195M/OC195M.methylation.summary.tsv" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/profile-fshd-methylation/OC195M/OC195M.methylation/4qA_all.stats.tsv" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/profile-fshd-methylation/OC195M/OC195M.methylation/4qA_all.methyl.bed.gz" .

aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/extract-classified-read-subsets/OC195M/OC195M.classified-subsets/OC195M.classified-subsets.manifest.tsv" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/extract-classified-read-subsets/OC195M/OC195M.classified-subsets/OC195M.4qA_all.bam" .
aws s3 cp "s3://alamyasingapore-nus-lab-production-processing/deployment-testing/OC195M_FSHD/extract-classified-read-subsets/OC195M/OC195M.classified-subsets/OC195M.4qA_all.bam.bai" .

