import "impute_sumstats_sub.wdl" as sub

task lift {

    String docker
    File sumstat_file
    File tbi_file = sumstat_file + ".tbi"
    String base = basename(sumstat_file)
    File b37_ref
    File b38_ref
    String dollar = "$"

    command <<<

        echo "COVID-19 HGI meta-analysis - lift over sumstats to b37"
        echo "${sumstat_file}"
        echo "${b37_ref}"
        echo "${b38_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${tbi_file} ${base}.tbi

        tabix -R ${b37_ref} ${base} | wc -l > b37.txt
        tabix -R ${b38_ref} ${base} | wc -l > b38.txt

        echo "`date` `cat b37.txt` chr 21 variants build 37"
        echo "`date` `cat b38.txt` chr 21 variants build 38"

        if ((`cat b37.txt` == 0 && `cat b38.txt` == 0)); then
            echo "`date` no chr 21 variants found in either build, quitting"
            exit 1
        fi

        if ((`cat b38.txt` > `cat b37.txt`)); then
            echo "`date` lifting to build 37"
            time /META_ANALYSIS/scripts/lift.py -chr "#CHR" -pos POS -ref Allele1 -alt Allele2 \
            -chain_file /liftover/hg38ToHg19.over.chain.gz -tmp_path /cromwell_root/ \
            ${base} > ${base}.lift.out 2> ${base}.lift.err
            gunzip -c ${base}.lifted.gz | \
            cut -f2- | awk '
            BEGIN { FS=OFS="\t" }
            NR==1 { for (i=1;i<=NF;i++) a[$i]=i; print $0 }
            NR>1 {
                temp=$a["#CHR"]; $a["#CHR"]=$a["anew_chr"]; $a["anew_chr"]=temp; temp=$a["POS"]; $a["POS"]=$a["anew_pos"]; $a["anew_pos"]=temp;
                sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]);
                if ($a["#CHR"] ~ /^[0-9]+$/) {
                    print $0
                }
            }' | bgzip > ${base}
        else
            echo "`date` presumably already in build 37"
        fi

    >>>

    output {
        File out = base
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

task harmonize {

    String docker
    File sumstat_file
    String base = basename(sumstat_file)
    File gnomad_ref
    String gnomad_ref_base = basename(gnomad_ref)
    Int n
    String options

    command <<<

        echo "COVID-19 HGI meta-analysis - harmonize sumstats to reference"
        echo "${sumstat_file}"
        echo "${gnomad_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${gnomad_ref} ${gnomad_ref_base}

        echo "`date` harmonizing stats with gnomAD"
        python3 /META_ANALYSIS/scripts/harmonize.py ${base} ${gnomad_ref_base} ${n} ${options}\
        | bgzip > ${base}.${gnomad_ref_base} && \
        tabix -S 1 -s 1 -b 2 -e 2 ${base}.${gnomad_ref_base} && \
        echo "`date` done"

    >>>

    output {
        File out = base + "." + gnomad_ref_base
        File out_tbi = base + "." + gnomad_ref_base + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}


task jass_preproc {

    String docker
    String pop
    String ref_panel_pop
    Float minMAF
    Float minAF
    Float maxAF
    File ref_panel_bim = sub(ref_panel_pop,"POP",pop)
    File sumstat_file
    String fname = basename(sumstat_file)
    String name
    
        
    command <<<
     
        echo "`date` filtering AF"

        echo "${sumstat_file}"
        echo "${fname}"

        mv ${sumstat_file} ${fname}
     
        gunzip -c ${fname} | awk '
        NR==1 {for (i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["AF_Allele2"] > ${minAF} && $a["AF_Allele2"] < ${maxAF}
        ' | gzip > ${fname}.AF_filtered.txt.gz

        ls -lh

        echo "COVID-19 HGI meta-analysis - jass preprocessing"

        python3 <<EOF
        import pandas as pd
        from collections import OrderedDict
        print("Writing GWAS info...")
        filename = '${fname}.AF_filtered.txt.gz'
        print(filename)
        d = OrderedDict([('filename', filename),
        ('internalDataLink', ''),
        ('Consortium', 'COVID'),
        ('Outcome', '${name}'),
        ('FullName', ''),
        ('Type', 'infection'),
        ('Reference', ''),
        ('ReferenceLink', ''),
        ('dataLink', ''),
        ('Nsample', int(filename.split('.')[7])+int(filename.split('.')[8])),
        ('Ncase', int(filename.split('.')[7])),
        ('Ncontrol', int(filename.split('.')[8])),
        ('Nsnp', ''),
        ('snpid', ''),
        ('a1', 'Allele1'),
        ('a2', 'Allele2'),
        ('freq', 'AF_Allele2'),
        ('pval', 'p.value'),
        ('n', 'N'),
        ('z', 'BETA'),
        ('OR', ''),
        ('se', 'SE'),
        ('CHR', '#CHR'),
        ('POS', 'POS'),
        ('code', ''),
        ('imp', ''),
        ('ncas', ''),
        ('ncont', ''),
        ('altNcas', ''),
        ('altNcont', ''),
        ('index_type', 'positional')])
        pd.DataFrame(data=d, index=[0]).to_csv('/cromwell_root/gwas_info.csv', sep='\t', index=False)
        EOF

        mkdir z_scores_chr
        
        jass_preprocessing --gwas-info gwas_info.csv \
        --ref-path ${ref_panel_bim} \
        --input-folder /cromwell_root/ \
        --diagnostic-folder /cromwell_root/ \
        --output-folder-1-file /cromwell_root/ \
        --output-folder z_scores_chr/ \
        --minimum-MAF ${minMAF}
      
        echo "`date` done"

    >>>

    output {
        File out_scores = "/cromwell_root/z_COVID_" + name + ".txt"
        Array[File] out_scores_chr = glob("/cromwell_root/z_scores_chr/z_*")
    }

    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 400 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow impute_sumstats {

    File sumstats_loc
    Array[Array[String]] sumstat_files = read_tsv(sumstats_loc)
    String gnomad_ref_template
    Array[String] chrom_list

    scatter (sumstat_file in sumstat_files) {

        call lift {
            input: sumstat_file=sumstat_file[0]
        }
        
        call harmonize {
            input: sumstat_file=lift.out, gnomad_ref=sub(gnomad_ref_template, "POP", sumstat_file[1]), n=sumstat_file[2]
        }

        call jass_preproc {
            input:
            sumstat_file=harmonize.out,
            pop=sumstat_file[1],
            name=sumstat_file[3]
        }

        call sub.raiss_combine{
            input:
            chrom_list=chrom_list,
            z_scores=jass_preproc.out_scores_chr,
            sumstat_file=sumstat_file,
            sumstat_lift=harmonize.out
        }

    }

}
