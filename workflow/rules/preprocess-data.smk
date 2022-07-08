rule get_functional_annotation_maf:
    input: "resources/PCAWG.public.maf.gz"
    output: "results/PCAWG.annotated.csv"
    shell:
        """
        gzip -dc {input} | 
            awk -F'\t' -vOFS=',' \
                'BEGIN{print "Gene", "Consequence","Substitution","Context","Tumor Type","Donor"}
                $7=="SNP" {print $1, $6, $8">"$10, substr($16,10,3), $42, $43}' \
                > {output}
        """