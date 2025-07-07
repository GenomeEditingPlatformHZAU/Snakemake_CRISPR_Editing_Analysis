rule HiTom_single_1_process_data_gunzip:
    input:
        "Resources/HiTom/{}".format(config["DATA"]["HiTom"]["reads"])
    output:
        temp(directory("{prefix}-Split"))
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"]
    shell:
        '''
        tar -xzvf {input}
        '''

rule HiTom_single_2_process_data_sample:
    input:
        "{prefix}-Split"
    output:
        fq1 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R1.fastq.gz",
        fq2 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R2.fastq.gz"
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{prefix}-{HiTom_sample}",
        outdir = lambda w, output: os.path.dirname(output.fq1),
    log:
        "logs/{prefix}/HiTom_single_process_data_sample_{prefix}-{HiTom_sample}.log"
    shell:
        '''
        gzip {input}/{params.sampleid}_R1.fastq > {log} 2>&1
        mv {input}/{params.sampleid}_R1.fastq.gz {output.fq1}

        gzip {input}/{params.sampleid}_R2.fastq > {log} 2>&1
        mv {input}/{params.sampleid}_R2.fastq.gz {output.fq2}
        '''

rule HiTom_single_3_prepare_description:
    input:
        "Resources/HiTom/{}".format(config["DATA"]["HiTom"]["description"])
    output:
        sg1 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA1_description.txt",
        sg2 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA2_description.txt"
    params:
        sampleid = "{HiTom_sample}"
    shell:
        '''
        cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print "-a "$2" -g "$3}}' > {output.sg1}
        cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print "-a "$4" -g "$5}}' > {output.sg2}
        '''

if config["module"]["HiTom"]["CRISPR"]:
    rule HiTom_single_4_CRISPResso2_CRISPR:
        input:
            fq1 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R1.fastq.gz",
            fq2 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R2.fastq.gz",
            sg1 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA1_description.txt",
            sg2 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA2_description.txt"
        output:
            sh = "results/{prefix}_HiTom_single/{prefix}_run_command/{HiTom_sample}.sh",
            CRISPResso1 = report(directory("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R1"), patterns=["{prefix}-{HiTom_sample}_9.Alleles_frequency_table_around_sgRNA_{name}.png", "{prefix}-{HiTom_sample}_3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Hi-Tom_CRISPR.rst", category="2-1. HiTom CRISPR"),
            #CRISPResso_ids1 = report("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R1/{prefix}-{HiTom_sample}_3b.Insertion_deletion_substitutions_size_hist.png", caption="../report/Hi-Tom_CRISPR_IDS.rst", category="2-1. HiTom CRISPR"),
            CRISPResso2 = report(directory("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R2"), patterns=["{prefix}-{HiTom_sample}_9.Alleles_frequency_table_around_sgRNA_{name}.png", "{prefix}-{HiTom_sample}_3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Hi-Tom_CRISPR.rst", category="2-1. HiTom CRISPR"),
            #CRISPResso_ids2 = report("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R2/{prefix}-{HiTom_sample}_3b.Insertion_deletion_substitutions_size_hist.png", caption="../report/Hi-Tom_CRISPR_IDS.rst", category="2-1. HiTom CRISPR")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{prefix}-{HiTom_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.CRISPResso1).split("/")[0:3]),
            extr = config["params"]["CRISPResso"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/HiTom_single_CRISPResso2_CRISPR_{prefix}-{HiTom_sample}.log"
        shell:
            '''
            set +e
            echo CRISPResso -r1 {params.workdir}/{input.fq1} --output_folder {params.workdir}/{params.outdir} {params.extr} \
             | cat - {input.sg1} | sed 'N;s/\\n/ /' > {output.sh}
            echo CRISPResso -r1 {params.workdir}/{input.fq2} --output_folder {params.workdir}/{params.outdir} {params.extr} \
             | cat - {input.sg2} | sed 'N;s/\\n/ /' >> {output.sh}
            bash {output.sh} > {params.workdir}/{log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R1/{params.sampleid}_9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R1/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R2/{params.sampleid}_9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R2/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.CRISPResso1} {output.CRISPResso2}" >> {params.workdir}/{log} 2>&1
                for i in `ls {params.outdir}/CRISPResso_on_{params.sampleid}-R1/9.Alleles_frequency_table_around* {params.outdir}/CRISPResso_on_{params.sampleid}-R2/9.Alleles_frequency_table_around* | xargs`
                do
                dir=`dirname ${{i}}`
                name=`basename ${{i}}`
                mv ${{i}} ${{dir}}/{params.sampleid}_${{name}}
                done
                mv {params.outdir}/CRISPResso_on_{params.sampleid}-R1/3b.Insertion_deletion_substitutions_size_hist.png {params.outdir}/CRISPResso_on_{params.sampleid}-R1/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.png
                mv {params.outdir}/CRISPResso_on_{params.sampleid}-R1/3b.Insertion_deletion_substitutions_size_hist.pdf {params.outdir}/CRISPResso_on_{params.sampleid}-R1/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.pdf
                mv {params.outdir}/CRISPResso_on_{params.sampleid}-R2/3b.Insertion_deletion_substitutions_size_hist.png {params.outdir}/CRISPResso_on_{params.sampleid}-R2/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.png
                mv {params.outdir}/CRISPResso_on_{params.sampleid}-R2/3b.Insertion_deletion_substitutions_size_hist.pdf {params.outdir}/CRISPResso_on_{params.sampleid}-R2/{params.sampleid}_3b.Insertion_deletion_substitutions_size_hist.pdf
                exit $exitcode
            fi
            '''

    rule HiTom_single_5_CRISPResso2_CRISPR_list:
        input:
            expand("results/{prefix}_HiTom_single/{prefix}_run_command/{HiTom_sample}.sh", prefix = config["prefix"], HiTom_sample = HiTom_samples)
        output:
            "results/{prefix}_HiTom_single/All_{prefix}_CRISPR_samples_HiTom.txt"
        log:
            "logs/{prefix}/HiTom_single_CRISPResso2_{prefix}_CRISPR.log"
        shell:
            '''
            cat {input} > {output}
            '''

if config["module"]["HiTom"]["ABE"]:
    rule HiTom_single_6_CRISPResso2_ABE:
        input:
            fq1 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R1.fastq.gz",
            fq2 = "results/{prefix}_HiTom_single/{prefix}-Split/{prefix}-{HiTom_sample}-R2.fastq.gz",
            sg1 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA1_description.txt",
            sg2 = "results/{prefix}_HiTom_single/{prefix}_description/{HiTom_sample}_sgRNA2_description.txt"
        output:
            sh = "results/{prefix}_HiTom_single/{prefix}_run_command/{HiTom_sample}.sh",
            CRISPResso1 = report(directory("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R1"), patterns=["{prefix}-{HiTom_sample}_9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Hi-Tom_ABE.rst", category="1. HiTom_ABE"),
            CRISPResso2 = report(directory("results/{prefix}_HiTom_single/{prefix}_results/CRISPResso_on_{prefix}-{HiTom_sample}-R2"), patterns=["{prefix}-{HiTom_sample}_9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Hi-Tom_ABE.rst", category="1. HiTom_ABE")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{prefix}-{HiTom_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.CRISPResso1).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/HiTom_single_CRISPResso2_ABE_{prefix}-{HiTom_sample}.log"
        shell:
            '''
            set +e
            echo CRISPResso -r1 {params.workdir}/{input.fq1} --output_folder {params.workdir}/{params.outdir} {params.extr} --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G \
             | cat - {input.sg1} | sed 'N;s/\\n/ /' > {output.sh}
            echo CRISPResso -r1 {params.workdir}/{input.fq2} --output_folder {params.workdir}/{params.outdir} {params.extr} --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G \
             | cat - {input.sg2} | sed 'N;s/\\n/ /' >> {output.sh}
            bash {output.sh} >> {params.workdir}/{log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R1/{params.sampleid}_9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {params.outdir}/CRISPResso_on_{params.sampleid}-R2/{params.sampleid}_9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.CRISPResso1} {output.CRISPResso2}" >> {params.workdir}/{log} 2>&1
                for i in `ls {params.outdir}/CRISPResso_on_{params.sampleid}-R1_{params.sampleid}-R2/9.Alleles_frequency_table_around* | xargs`
                do
                dir=`dirname ${{i}}`
                name=`basename ${{i}}`
                mv ${{i}} ${{dir}}/{params.sampleid}_${{name}}
                done
                exit $exitcode
            fi
            '''

    rule HiTom_single_7_CRISPResso2_ABE_list:
        input:
            expand("results/{prefix}_HiTom_single/{prefix}_run_command/{HiTom_sample}.sh", prefix = config["prefix"], HiTom_sample = HiTom_samples)
        output:
            "results/{prefix}_HiTom_single/All_{prefix}_ABE_samples_HiTom.txt"
        log:
            "logs/{prefix}/HiTom_single_CRISPResso2_{prefix}_ABE.log"
        shell:
            '''
            cat {input} > {output}
            '''
