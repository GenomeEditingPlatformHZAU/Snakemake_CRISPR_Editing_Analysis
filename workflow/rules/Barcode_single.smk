def get_read_file_path(wildcards):
  """ fetch the paths to the long read files from the configuration """
  return "Resources/Barcode/"+config["DATA"]["Barcode"]["fqs"][0], "Resources/Barcode/"+config["DATA"]["Barcode"]["fqs"][1]

rule Barcode_single_1_merge_fq:
    input:
        get_read_file_path
        #"Resources/Barcode/{}".format(config["DATA"]["Barcode"]["fqs"])
    output:
        "results/{prefix}_Barcode_single_Results/Split_fastqs/{prefix}.clean.fastq"
    threads: 1
    shell:
        '''
        zcat {input} > {output}
        '''

rule Barcode_single_2_process_barcode_primer:
    input:
        barcode = "results/Barcode_Primer/ER_{prefix}_barcode.txt" if config["module"]["Barcode_Design"] else "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["BarcodePrimer"])
    output:
        Fa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-F.fa"),
        FRa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-FR.fa"),
        FCa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-FC.fa"),
        FRCa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-FRC.fa"),
        FRP = "results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-F.txt",
        Ra = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-R.fa"),
        RRa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-RR.fa"),
        RCa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-RC.fa"),
        RRCa = temp("results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-RRC.fa"),
        RRP = "results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-R.txt"
    threads: 1
    conda:
        "../envs/Barcode_Hi-Tom.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/{prefix}/Barcode_single_process_barcode_primer_{prefix}_{Barcode_sample}.log"
    shell:
        '''
        cat {params.workdir}/{input.barcode} | sed 'N;s/\\n/\\t/' | sed "s/-F//g" | awk -F"\\t" '$1 == "{params.sampleid}" {{print ">"$1"\\n"$2}}' > {output.Fa}
        cat {params.workdir}/{input.barcode} | sed 'N;s/\\n/\\t/' | sed "s/-R//g" | awk -F"\\t" '$3 == "{params.sampleid}" {{print ">"$3"\\n"$4}}' > {output.Ra}

        seqkit seq -t DNA -r -w 200 {output.Fa} > {output.FRa}
        seqkit seq -t DNA -p -w 200 {output.Fa} > {output.FCa}
        seqkit seq -t DNA -r -p -w 200 {output.Fa} > {output.FRCa}
        seqkit seq -t DNA -r -w 200 {output.Ra} > {output.RRa}
        seqkit seq -t DNA -p -w 200 {output.Ra} > {output.RCa}
        seqkit seq -t DNA -r -p -w 200 {output.Ra} > {output.RRCa}
        cat {output.Fa} {output.FRa} {output.FCa} {output.FRCa} | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"\\t":$0 }}' | cut -f 2 | xargs | sed "s/ /\\n/g" > {output.FRP}
        cat {output.Ra} {output.RRa} {output.RCa} {output.RRCa} | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"\\t":$0 }}' | cut -f 2 | xargs | sed "s/ /\\n/g" > {output.RRP}
        '''

rule Barcode_single_3_process_barcode_primer_fq:
    input:
        FRP = "results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-F.txt",
        RRP = "results/{prefix}_Barcode_single_Results/Primers/{Barcode_sample}-R.txt",
        fqs = "results/{prefix}_Barcode_single_Results/Split_fastqs/{prefix}.clean.fastq"
    output:
        gzf = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-F.fastq.gz",
        gzr = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-R.fastq.gz"
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output.gzf)
    shell:
        '''
        set +e
        gzf_name=`basename {output.gzf}`
        gzr_name=`basename {output.gzr}`
        cd {params.outdir}

        cat {params.workdir}/{input.fqs} | \\grep -B1 -A2 -i -f {params.workdir}/{input.FRP} | grep -v "^--$" >> ${{gzf_name%.*}}
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            touch ${{gzf_name%.*}}
            gzip ${{gzf_name%.*}}
            exit 0
        else
            echo "Success run"
            gzip ${{gzf_name%.*}}

        fi
        ##################
        cat {params.workdir}/{input.fqs} | \\grep -B1 -A2 -i -f {params.workdir}/{input.RRP} | grep -v "^--$" >> ${{gzr_name%.*}}
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            touch ${{gzr_name%.*}}
            gzip ${{gzr_name%.*}}
            exit 0
        else
            echo "Success run"
            gzip ${{gzr_name%.*}}

        fi
        '''

rule Barcode_single_4_prepare_sample_description:
    input:
        "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["description"])
    output:
        Fd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-F_description.txt",
        Rd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-R_description.txt"
    params:
        workdir = config["workdir"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print "-a "$2" -g "$3}}' > {output.Fd}
        cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print "-a "$4" -g "$5}}' > {output.Rd}
        '''

if config["module"]["Barcode"]["CRISPR"]:
    rule Barcode_single_5_CRISPResso2_CRISPR_command:
        input:
            gzf = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-F.fastq.gz",
            gzr = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-R.fastq.gz",
            Fd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-F_description.txt",
            Rd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-R_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            echo CRISPResso -r1 {params.workdir}/{input.gzf} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/Cas9_Forward {params.extr} \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}

            echo CRISPResso -r1 {params.workdir}/{input.gzr} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/Cas9_Reverse {params.extr} \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            '''

    rule Barcode_single_6_CRISPResso2_CRISPR:
        input:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR.sh"
        output:
            FCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/Cas9_Forward/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR"),
            RCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/Cas9_Reverse/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_single_CRISPResso2_CRISPR_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {output.FCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.FCRISPResso}/3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.FCRISPResso} {output.RCRISPResso}" >> {params.workdir}/{log} 2>&1
                exit $exitcode
            fi
            '''

    rule Barcode_single_7_CRISPResso2_CRISPR_list:
        input:
            expand("results/{prefix}_Barcode_single_Results/Cas9_Forward/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples),
            expand("results/{prefix}_Barcode_single_Results/Cas9_Reverse/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples)
        output:
            "{prefix}_Barcode_single_Results/All_CRISPR_samples_Barcode.txt"
        shell:
            '''
            ls {input} > {output}
            '''

if config["module"]["Barcode"]["ABE"]:
    rule Barcode_single_8_CRISPResso2_ABE_command:
        input:
            gzf = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-F.fastq.gz",
            gzr = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-R.fastq.gz",
            Fd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-F_description.txt",
            Rd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-R_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            echo CRISPResso -r1 {params.workdir}/{input.gzf} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/ABE_Forward --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            
            echo CRISPResso -r1 {params.workdir}/{input.gzr} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/ABE_Reverse --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} \
             | cat - {input.Rd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            '''

    rule Barcode_single_9_CRISPResso2_ABE:
        input:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE.sh"
        output:
            FCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/ABE_Forward/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE"),
            RCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/ABE_Reverse/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_single_CRISPResso2_ABE_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {output.FCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.FCRISPResso} {output.RCRISPResso}" >> {params.workdir}/{log} 2>&1
                exit $exitcode
            fi
            '''

    rule Barcode_single_10_CRISPResso2_ABE_list:
        input:
            expand("results/{prefix}_Barcode_single_Results/ABE_Forward/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples),
            expand("results/{prefix}_Barcode_single_Results/ABE_Reverse/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples)
        output:
            "results/{prefix}_Barcode_single_Results/All_ABE_samples_Barcode.txt"
        shell:
            '''
            ls {input} > {output}
            '''

if config["module"]["Barcode"]["CBE"]:
    rule Barcode_single_11_CRISPResso2_CBE_command:
        input:
            gzf = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-F.fastq.gz",
            gzr = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-R.fastq.gz",
            Fd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-F_description.txt",
            Rd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-R_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            echo CRISPResso -r1 {params.workdir}/{input.gzf} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/CBE_Forward --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            
            echo CRISPResso -r1 {params.workdir}/{input.gzr} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/CBE_Reverse --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} \
             | cat - {input.Rd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            '''

    rule Barcode_single_12_CRISPResso2_CBE:
        input:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE.sh"
        output:
            FCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/CBE_Forward/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE"),
            RCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/CBE_Reverse/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_single_CRISPResso2_CBE_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {output.FCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.FCRISPResso} {output.RCRISPResso}" >> {params.workdir}/{log} 2>&1
                exit $exitcode
            fi
            '''

    rule Barcode_single_13_CRISPResso2_CBE_list:
        input:
            expand("results/{prefix}_Barcode_single_Results/CBE_Forward/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples),
            expand("results/{prefix}_Barcode_single_Results/CBE_Reverse/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples)
        output:
            "results/{prefix}_Barcode_single_Results/All_CBE_samples_Barcode.txt"
        shell:
            '''
            ls {input} > {output}
            '''

if config["module"]["Barcode"]["cpf1"]:
    rule Barcode_single_14_CRISPResso2_cpf1_command:
        input:
            gzf = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-F.fastq.gz",
            gzr = "results/{prefix}_Barcode_single_Results/Split_fastqs/{Barcode_sample}-R.fastq.gz",
            Fd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-F_description.txt",
            Rd = "results/{prefix}_Barcode_single_Results/description/{Barcode_sample}-R_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            echo CRISPResso -r1 {params.workdir}/{input.gzf} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/cpf1_Forward {params.extr} --cleavage_offset -6 --plot_window_size 25 --default_min_aln_score 50 \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}

            echo CRISPResso -r1 {params.workdir}/{input.gzr} --name {params.sampleid} --output_folder {params.workdir}/{params.outdir}/cpf1_Reverse {params.extr} --cleavage_offset -6 --plot_window_size 25 --default_min_aln_score 50 \
             | cat - {input.Fd} | sed 'N;s/\\n/ /' >> {output.bash_command}
            '''

    rule Barcode_single_15_CRISPResso2_cpf1:
        input:
            bash_command = "results/{prefix}_Barcode_single_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1.sh"
        output:
            FCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/cpf1_Forward/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_cpf1.rst", category="2-1. Barcode cpf1"),
            RCRISPResso = report(directory("results/{prefix}_Barcode_single_Results/cpf1_Reverse/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_cpf1.rst", category="2-1. Barcode cpf1")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_single_CRISPResso2_cpf1_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {output.FCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.FCRISPResso}/3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.RCRISPResso}/3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                exit 0
            else
                echo "Success run {output.FCRISPResso} {output.RCRISPResso}" >> {params.workdir}/{log} 2>&1
                exit $exitcode
            fi
            '''

    rule Barcode_single_16_CRISPResso2_cpf1_list:
        input:
            expand("results/{prefix}_Barcode_single_Results/cpf1_Forward/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples),
            expand("results/{prefix}_Barcode_single_Results/cpf1_Reverse/CRISPResso_on_{Barcode_sample}", prefix = config["prefix"], Barcode_sample = Barcode_samples)
        output:
            "{prefix}_Barcode_single_Results/All_cpf1_samples_Barcode.txt"
        shell:
            '''
            ls {input} > {output}
            '''
