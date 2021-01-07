task extract_af {

    File vcf
    String base_vcf = basename(vcf)

    command <<<

        mv ${vcf} ${base_vcf}

        python3 - ${base_vcf} <<EOF

        import sys
        import gzip

        #v3
        #POPS = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas', 'oth']
        POPS = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'nfe_est', 'nfe_nwe', 'nfe_seu', 'nfe_onf', 'oth']

        out_tsv = [gzip.open('ref_${base_vcf}_' + pop + '.gz', 'wt') for pop in POPS]
        for out in out_tsv:
            out.write('#chr\tpos\tref\talt\taf_alt\tfilter\tan\n')

        af_fields = ['AF_' + pop for pop in POPS]
        with gzip.open(sys.argv[1], 'rt') as f:
            line = f.readline().strip()
            while line.startswith('##'):
                line = f.readline().strip()
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chrom = s[0].replace('chr', '').replace('X', '23')
                cpra = chrom + '\t' + s[1] + '\t' + s[3] + '\t' + s[4]
                info = {dat.split('=')[0]: dat.split('=')[1] for dat in s[7].split(';') if '=' in dat}
                for i,out in enumerate(out_tsv):
                    af = float(info[af_fields[i]]) if af_fields[i] in info else 'NA'
                    out.write(cpra + '\t' + str(af) + '\t' + s[6] + '\t' + info['AN'] + '\n')

        for out in out_tsv:
            out.close()

        EOF

    >>>

    output {
        Array[File] out = glob('ref_*')
    }

    runtime {
        docker: "gcr.io/covid-19-hg/meta:1d50c"
        cpu: 1
        # 40M variants in variant_file to look up takes about 4G
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task combine {

    Array[Array[File]] files2D
    Array[File] files = flatten(files2D)
    #v3
    #Array[String] pops = ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "sas", "oth"]
    Array[String] pops = ["afr", "amr", "asj", "eas", "fin", "nfe", "nfe_est", "nfe_nwe", "nfe_seu", "nfe_onf", "oth"]
    String dollar = "$"

    command <<<

        for file in ${sep=" " files}; do
            mv $file `basename $file`
        done
        while read pop; do
            files=(*_$pop*)
            cat <(gunzip -c ${dollar}{files[0]} | head -1) \
                <(for file in ${dollar}{files[@]}; do gunzip -c $file | tail -n+2; done) \
            | bgzip > gnomad_v2.1.1_b37_ref_$pop.gz
            tabix -s 1 -b 2 -e 2 gnomad_v2.1.1_b37_ref_$pop.gz
        done < ${write_lines(pops)}

    >>>

    output {
        Array[File] out = glob("*.gz")
        Array[File] tbi = glob("*.tbi")
    }

    runtime {
        docker: "gcr.io/covid-19-hg/meta:1d50c"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 200 SSD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow gnomad_extract_af {

    File vcf_files
    Array[String] vcfs = read_lines(vcf_files)

    scatter (vcf in vcfs) {
        call extract_af {
            input: vcf=vcf
        }
    }

    call combine {
        input: files2D=extract_af.out
    }
}
