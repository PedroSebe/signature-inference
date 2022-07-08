import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

rule download_everything:
    input:
        "resources/genome/hg19.fa",
        "resources/genome/hg19.fa.fai",
        "resources/COSMIC_v3.2_SBS_GRCh37.txt",
        "resources/pancan_pcawg_2020/data_mutations.txt",
        "resources/msk_impact_2017/data_mutations.txt",
        "resources/Target_TSO.bed",
        "resources/Target_MSK.bed",
        "resources/PCAWG.public.maf.gz",
        "results/NikZainal-breast_cancer_spectra.csv"

rule reference_genome:
    input: HTTP.remote("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz")
    output: "resources/genome/hg19.fa"
    threads: 1
    shell:
        """
        mkdir -p resources/genome
        mv {input} {output}.gz
        gunzip {output}.gz
        """
rule reference_genome_index:
    input: "resources/genome/hg19.fa"
    output: "resources/genome/hg19.fa.fai"
    threads: 1
    shell: "samtools faidx {input}"

rule cosmic_signatures:
    input: HTTP.remote("https://cancer.sanger.ac.uk/signatures/documents/452/COSMIC_v3.2_SBS_GRCh37.txt")
    output: "resources/COSMIC_v3.2_SBS_GRCh37.txt"
    shell: "mv {input} {output}"

rule pcawg_all_data:
    input: HTTP.remote("https://cbioportal-datahub.s3.amazonaws.com/pancan_pcawg_2020.tar.gz")
    output: "resources/pancan_pcawg_2020/data_mutations.txt"
    shell: "tar -xzf {input} -C resources/"

rule msk_all_data:
    input: HTTP.remote("https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz")
    output: "resources/msk_impact_2017/data_mutations.txt"
    shell: "tar -xzf {input} -C resources/"

rule TSO_target_BED:
    input: HTTP.remote("https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/nextera/nextera-flex-for-enrichment/TruSight_One_TargetedRegions_v1.1.bed")
    output: "resources/Target_TSO.bed"
    shell: "mv {input} {output}"

rule MSK_target_spreadsheet:
    input: HTTP.remote("https://www.jmdjournal.org/cms/10.1016/j.jmoldx.2014.12.006/attachment/b2f0836d-3bbb-4fff-ab49-8765fb647a73/mmc4.xlsx")
    output: "resources/Target_MSK.xlsx"
    shell: "mv {input} {output}"

rule MSK_target_BED:
    input:
        excel = "resources/Target_MSK.xlsx",
        script = "workflow/scripts/generate_msk_bed.py"
    output: "resources/Target_MSK.bed"
    shell: "python {input.script} {input.excel} > {output}"

rule pcawg_spectra:
    input: HTTP.remote("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Input_Data_PCAWG7_23K_Spectra_DB/Mutation_Catalogs_--_Spectra_of_Individual_Tumours/WGS_PCAWG_2018_02_09.zip")
    output:
        "resources/PCAWG_catalogues/WGS_PCAWG.1536.csv",
        "resources/PCAWG_catalogues/WGS_PCAWG.96.csv",
        "resources/PCAWG_catalogues/WGS_PCAWG.indels.csv",
        "resources/PCAWG_catalogues/WGS_PCAWG.192.csv",
        "resources/PCAWG_catalogues/WGS_PCAWG.dinucs.csv"
    shell:
        """
        mkdir -p resources/PCAWG_catalogues
        unzip {input} -d resources/PCAWG_catalogues
        rm -r resources/PCAWG_catalogues/__MACOSX
        """

rule pcawg_maf:
    input: HTTP.remote("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz")
    output: "resources/PCAWG.public.maf.gz"
    shell: "mv {input} {output}"

rule NikZainal_FTP:
    input: FTP.remote("ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl/SUBSTITUTIONS_13Apr2012_snz.txt")
    output: "resources/NikZainal-breast_cancer.csv"
    shell: "mv {input} {output}"

rule NikZainal_spectra:
    input: "resources/NikZainal-breast_cancer.csv"
    output: "results/NikZainal-breast_cancer_spectra.csv"
    shell:
        """
        awk 'NR>1{{print $2, substr($8, 10, 1)"["$9">"$10"]"substr($11, 1, 1)}}' {input} |
            sort | uniq -c | 
            awk -vOFS='\t' \
                'BEGIN{{print "sample", "substitution_type", "count"}} {{print $2, $3, $1}}' \
                > {output}
        """