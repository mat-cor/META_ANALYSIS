task raiss {

    String docker
    String chrom
    String pop
    String name
    String gwas = "COVID_" + name
    
    Array[File] z_scores

    String ref_panel_chr_pop
    String ref_panel_chr=sub(ref_panel_chr_pop,"POP",pop)
    File ref=sub(ref_panel_chr,"CHROM",chrom)
    
    String ref_panel_path_pop
    String ref_panel_path=sub(ref_panel_path_pop,"POP",pop)
        
    
    String ld_slices_path_pop
    String ld_slices_path=sub(ld_slices_path_pop,"POP",pop)

    String ld_slices_chr_pop
    String ld_slices_chr=sub(ld_slices_chr_pop,"POP",pop)
    File ld_slices_file=sub(ld_slices_chr,"CHROM",chrom)
    Array[Array[File]] ld_slices=read_tsv(ld_slices_file)
    
    Float eigen_threshold
    Float r2_threshold

    command <<<

        echo "COVID-19 HGI meta-analysis - RAISS sumstats imputation"

        mkdir z_scores
        mkdir z_scores_imputed
        
        find . -name "z_${gwas}_${chrom}.txt" -exec cp {} ./z_scores/ \;
              
        raiss --chrom ${chrom} \
        --gwas ${gwas} \
        --ref-folder /cromwell_root/${ref_panel_path} \
        --ld-folder /cromwell_root/${ld_slices_path} \
        --zscore-folder z_scores \
        --output-folder z_scores_imputed \
        --ld-type scipy \
        --ref-panel-suffix .bim \
        --eigen-threshold ${eigen_threshold} \
        --R2-threshold ${r2_threshold}

        mkdir z_scores_masked
        mkdir z_scores_imputed_test
        mkdir imputation_performance

        echo 'null' > imputation_performance/null.csv

        if [[ ${chrom} == *21 ]] || [[ ${chrom} == *22 ]]
        then
        python3 <<EOF
        import raiss
        import pandas as pd
        base_dir = '/cromwell_root'
        code_GWAS = '${gwas}'
        chr = '${chrom}'
        perf_results = raiss.imputation_R2.grid_search("{0}/z_scores/".format(base_dir),
                                        "{0}/z_scores_masked/".format(base_dir),
                                        "{0}/z_scores_imputed_test/".format(base_dir),
                                        "{0}/${ref_panel_path}/".format(base_dir),
                                        "{0}/${ld_slices_path}/".format(base_dir),
                                        code_GWAS,
                                        eigen_ratio_grid =[0.1, 0.001, 0.000001,0.00000001], 
                                        N_to_mask=20000,
                                        ref_panel_suffix = ".bim",
                                        chrom=chr,
                                        ld_type="scipy")
        perf_results.to_csv("{0}/imputation_performance/performance_{1}_{2}.csv".format(base_dir, code_GWAS, chr))
        EOF
        fi

        echo "`date` done"

    >>>

    output {
        File out_imp_scores = "/cromwell_root/z_scores_imputed/z_" + gwas + "_" + chrom + ".txt"
        Array[File] out_test = glob("/cromwell_root/imputation_performance/*.csv")
    }

    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}


task combine {

    String docker
    File sumstat_lift
    String out_prefix
    Array[File] results
    
    command <<<

        cat \
        <(head -n 1 ${results[0]}) \
        <(for file in ${sep=" " results}; do tail -n+2 $file; done) \
        | bgzip > ${out_prefix}.imputed.txt.gz


        Rscript /META_ANALYSIS/scripts/combine_zscores.R \
        -s ${sumstat_lift} \
        -i ${out_prefix}.imputed.txt.gz \
        -o ${out_prefix}.combined_imputed.txt

        bgzip ${out_prefix}.combined_imputed.txt


    >>>

    output {

        File z_scores_imputed = "${out_prefix}.imputed.txt.gz"
        File out = "${out_prefix}.combined_imputed.txt.gz"
        
    }

    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "16 GB"
        disks: "local-disk 50 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}


task plot {

    File sumstat
    String base = basename(sumstat)
    String docker
    Int loglog_ylim

    command <<<

        mv /cromwell_root/covid19-hg-cromwell/impute_sumstats/*/call-raiss_combine/shard-*/sub.raiss_combine/*/call-combine/*.combined_imputed.txt.gz /cromwell_root/

        qqplot.raiss.R --file ${base} --bp_col "POS" --chrcol "#CHR" --pval_col "p.value" --snp_col "rsID" --loglog_ylim ${loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: 10*ceil(size(sumstat, "G")) + " GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}



workflow raiss_combine {

    Array[String] chrom_list
    File sumstat_lift
    Array[String] sumstat_file
    Array[File] z_scores

    scatter (chrom in chrom_list){ 
        call raiss {
            input:
            chrom = chrom,
            name = sumstat_file[3],
            z_scores = z_scores,
            pop=sumstat_file[1]
        } 
    }

    call combine {
        input:
        sumstat_lift = sumstat_lift,
        out_prefix = sumstat_file[3],
        results = raiss.out_imp_scores
    }

    call plot {
        input:
        sumstat = combine.out
    }

    output {
        File out = combine.out
    }
}