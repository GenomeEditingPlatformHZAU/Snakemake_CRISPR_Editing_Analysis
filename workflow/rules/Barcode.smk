def get_read_file_path(wildcards):
  """ fetch the paths to the long read files from the configuration """
  return "Resources/Barcode/"+config["DATA"]["Barcode"]["fqs"][0], "Resources/Barcode/"+config["DATA"]["Barcode"]["fqs"][1]

################################################################################################################################################################
########################################################## 1. 处理测序fq数据 ####################################################################################
################################################################################################################################################################
rule Barcode_1_flash_merge_fq:
    input:
        get_read_file_path
    output:
        "results/{prefix}_Barcode_pair_Results/Split_fastqs/out.extendedFrags.fastq"
    threads: 2
    conda:
        "../envs/Barcode_Hi-Tom.yaml"
    params:
        outdir = lambda w, output: os.path.dirname(output[0]),
        extr = config["params"]["flash"]
    shell:
        '''
        flash {params.extr} -t {threads} --output-directory {params.outdir} {input}
        '''

rule Barcode_2_cat_merge_fq_skip_flash:
    input:
        get_read_file_path
    output:
        "results/{prefix}_Barcode_pair_Results/Split_fastqs/{prefix}.clean.fastq"
    shell:
        '''
        zcat {input} > {output}
        '''

################################################################################################################################################################
########################################################## 2. 解析样品引物文件 ##################################################################################
################################################################################################################################################################
rule Barcode_3_process_barcode_primer:
    '''
    perl -nle'BEGIN {{
            @map{{ A, a, C, c, G, g, T, t }} = ( T, t, G, g, C, c, A, a )
            }}
            print /^>/ ?
                $_ :
                    join //, map $map{{ $_ }}, split //, scalar reverse
            '  {output.RP} | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"\\t":$0 }}' | sed "s/^>//g" > {output.RCP}
    '''
    input:
        barcode = "results/Barcode_Primer/ER_{prefix}_barcode.txt" if config["module"]["Barcode_Design"] else "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["BarcodePrimer"])
    output:
        FP = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-F.txt"),
        RP = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-R.txt"),
        RCP = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-RC.txt"),
        Fa = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-F.fa"),
        Ra = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-R.fa"),
        Ca = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-C.fa"),
        RCa = temp("results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}-RC.fa"),
        FRP = "results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}.txt"
    threads: 1
    conda:
        "../envs/Barcode_Hi-Tom.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output["FRP"])
    shell:
        '''
        cat {input.barcode} | sed 'N;s/\\n/\\t/' | sed "s/-F//g" | awk -F"\\t" '$1 == "{params.sampleid}" {{print $1"\\t"$2}}' > {output.FP}
        cat {input.barcode} | sed 'N;s/\\n/\\t/' | sed "s/-R//g" | awk -F"\\t" '$3 == "{params.sampleid}" {{print ">"$3"\\n"$4}}' > {output.RP}
        python workflow/scripts/reverse_complement_fasta.py -i {output.RP} -o {output.RCP}
        cat {output.FP} {output.RCP} | sed 'N;s/\\n/ /' | awk '{{print $1,$2,$4}}' | awk 'BEGIN{{OFS="\\n"}} $1 == "{params.sampleid}" {{print ">"$1,$2"Y"$3}}' > {output.Fa}
        seqkit seq -t DNA -r -w 200 {output.Fa} > {output.Ra}
        seqkit seq -t DNA -p -w 200 {output.Fa} > {output.Ca}
        seqkit seq -t DNA -r -p -w 200 {output.Fa} > {output.RCa}
        sed -i "s/Y/\\.\\*/g" {output.Fa}
        sed -i "s/Y/\\.\\*/g" {output.Ra}
        sed -i "s/R/\\.\\*/g" {output.Ca}
        sed -i "s/R/\\.\\*/g" {output.RCa}

        cat {output.Fa} {output.Ra} {output.Ca} {output.RCa} | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"\\t":$0 }}' | cut -f 2 | xargs | sed "s/ /\\n/g" > {output.FRP}
        '''

################################################################################################################################################################
########################################################## 3. 拆分样品fq文件 ####################################################################################
################################################################################################################################################################
# rule Barcode_4_process_barcode_primer_fq:
#     input:
#         barcode_primer = "results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}.txt",
#         fqs = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{prefix}.clean.fastq" if config["params"]["skip_flash"] else "results/{prefix}_Barcode_pair_Results/Split_fastqs/out.extendedFrags.fastq"
#     output:
#         gz = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz"
#     threads: 1
#     params:
#         workdir = config["workdir"],
#         prefix = config["prefix"],
#         sampleid = "{Barcode_sample}",
#         outdir = lambda w, output: os.path.dirname(output.gz)
#     shell:
#         '''
#         set +e
#         gz_name=`basename {output.gz}`
#         cd {params.outdir}

#         cat {params.workdir}/{input.fqs} | \\grep -B1 -A2 -i -f {params.workdir}/{input.barcode_primer} | grep -v "^--$" >> ${{gz_name%.*}}
#         exitcode=$?
#         if [ $exitcode -ne 0 ]
#         then
#             touch ${{gz_name%.*}}
#             gzip ${{gz_name%.*}}
#             exit 0
#         else
#             echo "Success run"
#             gzip ${{gz_name%.*}}
#             exit $exitcode
#         fi
#         '''

rule Barcode_4_process_barcode_primer_fq:
    """
    若扩增序列长度小于120bp则跳过R1和R2的拼接过程;
    """
    input:
        description = "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["description"]),
        barcode_primer = "results/{prefix}_Barcode_pair_Results/Primers/{Barcode_sample}.txt",
        fqs_skip = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{prefix}.clean.fastq",
        fqs = "results/{prefix}_Barcode_pair_Results/Split_fastqs/out.extendedFrags.fastq"
    output:
        gz = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz"
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output.gz)
    shell:
        '''
        set +e
        gz_name=`basename {output.gz}`
        cd {params.outdir}

        # 查找匹配行并处理
        while IFS=$'\\t' read -r col1 col2 col3; do
            if [[ "$col1" == "{params.sampleid}" ]]; then
                len=${{#col2}}
                echo "匹配样本：$col1"
                echo "扩增序列长度：$len"
                if [ "$len" -lt 120 ]; then
                    echo "using {input.fqs_skip} file"
                    cat {params.workdir}/{input.fqs_skip} | \\grep -B1 -A2 -i -f {params.workdir}/{input.barcode_primer} | grep -v "^--$" >> ${{gz_name%.*}}
                else
                    echo "using {input.fqs} file"
                    cat {params.workdir}/{input.fqs} | \\grep -B1 -A2 -i -f {params.workdir}/{input.barcode_primer} | grep -v "^--$" >> ${{gz_name%.*}}
                fi
            fi
        done < {params.workdir}/{input.description}

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            touch ${{gz_name%.*}}
            gzip ${{gz_name%.*}}
            exit 0
        else
            echo "Success run"
            gzip ${{gz_name%.*}}
            exit $exitcode
        fi
        '''

rule Barcode_5_process_barcode_primer_fq_subgenome:
    input:
        fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
        description = "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["description"])
    output:
        directory("results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq")
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        sampleid = "{Barcode_sample}"
    shell:
        '''
        mkdir -p {output}
        subgenome=-1
        mysubgenome=(A D)
        nacol=`cat {input.description} | \\grep -w "{params.sampleid}" | cut -f 4`
        if [[ "${{nacol}}" == "NA" ]]; then
            echo "Normal"
            cp {input.fq} {output}/{params.sampleid}_subgenomeA.fastq.gz
        else
            echo "Subgenome"
            for i in `cat {input.description} | \\grep -w "{params.sampleid}" | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $5,$7}}' | sed "s/\\t/\\n/g" | xargs`
            do
            subgenome=$((subgenome + 1))
            zcat {input.fq} | \\grep -B1 -A2 -i "${{i}}" | grep -v "^--$" > {output}/{params.sampleid}_subgenome${{mysubgenome[${{subgenome}}]}}.fastq
            gzip {output}/{params.sampleid}_subgenome${{mysubgenome[${{subgenome}}]}}.fastq
            done
        fi
        '''

################################################################################################################################################################
########################################################## 4. 解析扩增序列和sgRNA ###############################################################################
################################################################################################################################################################
rule Barcode_6_prepare_sample_description:
    """
    为便于后续统计每个样品的每个靶标的总编辑率, 此处若某个样有两个sgRNA时分开单独跑, 且将结果存放于同一个文件夹下, 可能会对部分结果的统计有影响, 
    但如果只关心图9则没有影响, 因为其文件名称附加了对应的sgRNA序列;
    """
    input:
        "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["description"])
    output:
        "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
    params:
        workdir = config["workdir"],
        sampleid = "{Barcode_sample}",
        differentiation_subgenomes = config["params"]["differentiation_subgenomes"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        #cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print "-a "$2" -g "$3}}' > {output}
        #cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} {{if($1 == a && $3 ~ /,/) {{split($3,b,","); print "-a "$2" -g "b[1]"\\n-a "$2" -g "b[2];}} else if ($1 == a && $3 !~ /,/) {{print "-a "$2" -g "$3}}}}' > {output}
        #cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print $2,$3}}' | sed "s/,/\\t/g" | awk -F"\\t" '{{for(i=2;i<=NF;i++)print "-a "$1" -g "$i}}' > {output}
        if [[ "{params.differentiation_subgenomes}" == "True" ]]; then
            echo "Distinguishing the AD subgenome"
            cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} {{if ($1==a && $4=="NA") {{print a"_subgenomeA",$2,$3}} else if ($1==a && $4!="NA") {{print a"_subgenomeA",$4,$3"\\n"a"_subgenomeD",$6,$3}}}}' | sed "s/,/\\t/g" | awk -F"\\t" '{{for(i=3;i<=NF;i++)print "--name "$1" -a "$2" -g "$i}}' > {output}
        else
            echo "Run normal module"
            cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} $1 == a {{print a,$2,$3}}' | sed "s/,/\\t/g" | awk -F"\\t" '{{for(i=3;i<=NF;i++)print "--name "$1" -a "$2" -g "$i}}' > {output}
        fi
        '''

rule Barcode_7_prepare_sample_description_PE:
    input:
        "Resources/Barcode/{}".format(config["DATA"]["Barcode"]["description"])
    output:
        "results/{prefix}_Barcode_pair_Results/description_PE/{Barcode_sample}_description.txt"
    params:
        workdir = config["workdir"],
        sampleid = "{Barcode_sample}",
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        cat {input} | awk -F"\\t" -v a="{params.sampleid}" 'BEGIN{{OFS="\\t"}} {{if($1 == a && NF == 5) {{print "-a "$2" --prime_editing_pegRNA_extension_seq "$3" --prime_editing_pegRNA_spacer_seq "$4" --prime_editing_nicking_guide_seq "$5;}} else if ($1 == a && NF == 4) {{print "-a "$2" --prime_editing_pegRNA_extension_seq "$3" --prime_editing_pegRNA_spacer_seq "$4}}}}' > {output}
        '''

################################################################################################################################################################
########################################################## 5. 样品编辑检测 ######################################################################################
################################################################################################################################################################
if config["module"]["Barcode"]["CRISPR"]:
    rule Barcode_8_CRISPResso2_CRISPR_command:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            command="CRISPResso -r1 {params.workdir}/{input.fq} --output_folder {params.workdir}/{params.outdir}/Cas9 {params.extr} "
            sed "s|^|${{command}} |" {input.description} > {output.bash_command}
            '''

    rule Barcode_9_CRISPResso2_CRISPR_command_subgenome:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR_subgenome.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            for i in `ls {input.fq} | sort -V | xargs`
            do
            subgenome=`echo ${{i}} | sed "s/\\.fastq\\.gz//"`
            command="CRISPResso -r1 {params.workdir}/{input.fq}/${{i}} --output_folder {params.workdir}/{params.outdir}/Cas9_subgenome {params.extr} "
            cat {input.description} | awk -v a="${{subgenome}}" -v command="${{command}}" '$2==a {{print command,$0}}' >> {output.bash_command}
            done
            '''

    rule Barcode_10_CRISPResso2_CRISPR:
        input:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR_subgenome.sh" if config["params"]["differentiation_subgenomes"] else "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CRISPR.sh",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            CRISPResso_png = report(directory("results/{prefix}_Barcode_pair_Results/Cas9_subgenome/CRISPResso_on_{Barcode_sample}_subgenomeA"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR figure") if config["params"]["differentiation_subgenomes"] else report(directory("results/{prefix}_Barcode_pair_Results/Cas9/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png", "3b.Insertion_deletion_substitutions_size_hist.png"], caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR figure"),
            CRISPResso = report("results/{{prefix}}_Barcode_pair_Results/CRISPR_Alleles_frequency_subgenome/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR table") if config["params"]["differentiation_subgenomes"] else report("results/{{prefix}}_Barcode_pair_Results/CRISPR_Alleles_frequency/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_CRISPR.rst", category="2-1. Barcode CRISPR table")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            zaiti = config["zaiti"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_CRISPResso2_CRISPR_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.workdir}/{output.CRISPResso_png}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {params.workdir}/{output.CRISPResso_png}/3b.Insertion_deletion_substitutions_size_hist.png >> {params.workdir}/{log} 2>&1
                touch {output.CRISPResso}
                exit 0
            else
                echo "Success run {output.CRISPResso}" >> {params.workdir}/{log} 2>&1
                rows=`cat {input.description} | cut -f 6 -d " " | wc -l`
                for (( i=1; i<=$rows; i++ ))
                do
                    sgRNA=`cat {input.description} | cut -f 6 -d " " | sed -n ${{i}}p`
                    cat {params.workdir}/{output.CRISPResso_png}/Alleles_frequency_table_around_sgRNA_${{sgRNA}}.txt | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $3 ~ /False/ && $8 >= 0.2 {{print $0}}' | awk -F"\\t" -v a="{params.sampleid}" -v b="{params.prefix}" -v c="{params.zaiti}" -v d="sgRNA${{i}}" '{{SUM+=$8}}END{{print b"_"a"\\t"SUM"\\t"c"\\t"d}}' >> {output.CRISPResso}
                done
                exit $exitcode
            fi
            '''

    rule Barcode_11_CRISPResso2_CRISPR_list:
        input:
            expand("results/{prefix}_Barcode_pair_Results/CRISPR_Alleles_frequency_subgenome/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"]) if config["params"]["differentiation_subgenomes"] else expand("results/{prefix}_Barcode_pair_Results/CRISPR_Alleles_frequency/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"])
        output:
            "results/{prefix}_Barcode_pair_Results/All_CRISPR_samples_Barcode.txt"
        shell:
            '''
            cat {input} > {output}
            '''

if config["module"]["Barcode"]["ABE"]:
    rule Barcode_12_CRISPResso2_ABE_command:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            command="CRISPResso -r1 {params.workdir}/{input.fq} --output_folder {params.workdir}/{params.outdir}/ABE --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} "
            sed "s|^|${{command}} |" {input.description} > {output.bash_command}
            '''

    rule Barcode_13_CRISPResso2_ABE_command_subgenome:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE_subgenome.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            for i in `ls {input.fq} | sort -V | xargs`
            do
            subgenome=`echo ${{i}} | sed "s/\\.fastq\\.gz//"`
            command="CRISPResso -r1 {params.workdir}/{input.fq}/${{i}} --output_folder {params.workdir}/{params.outdir}/ABE_subgenome --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G {params.extr} "
            cat {input.description} | awk -v a="${{subgenome}}" -v command="${{command}}" '$2==a {{print command,$0}}' >> {output.bash_command}
            done
            '''

    rule Barcode_14_CRISPResso2_ABE:
        input:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE_subgenome.sh" if config["params"]["differentiation_subgenomes"] else "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_ABE.sh",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            CRISPResso_png = report(directory("results/{prefix}_Barcode_pair_Results/ABE_subgenome/CRISPResso_on_{Barcode_sample}_subgenomeA"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE figure") if config["params"]["differentiation_subgenomes"] else report(directory("results/{prefix}_Barcode_pair_Results/ABE/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE figure"),
            CRISPResso = report("results/{{prefix}}_Barcode_pair_Results/ABE_Alleles_frequency_subgenome/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE table") if config["params"]["differentiation_subgenomes"] else report("results/{{prefix}}_Barcode_pair_Results/ABE_Alleles_frequency/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_ABE.rst", category="2-1. Barcode ABE table")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            zaiti = config["zaiti"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_CRISPResso2_ABE_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.workdir}/{output.CRISPResso_png}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.CRISPResso}
                exit 0
            else
                echo "Success run {output.CRISPResso_png}" >> {params.workdir}/{log} 2>&1
                rows=`cat {input.description} | cut -f 6 -d " " | wc -l`
                for (( i=1; i<=$rows; i++ ))
                do
                    sgRNA=`cat {input.description} | cut -f 6 -d " " | sed -n ${{i}}p`
                    cat {params.workdir}/{output.CRISPResso_png}/Alleles_frequency_table_around_sgRNA_${{sgRNA}}.txt | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $3 ~ /False/ && $8 >= 0.2 {{print $0}}' | awk -F"\\t" -v a="{params.sampleid}" -v b="{params.prefix}" -v c="{params.zaiti}" -v d="sgRNA${{i}}" '{{SUM+=$8}}END{{print b"_"a"\\t"SUM"\\t"c"\\t"d}}' >> {output.CRISPResso}
                done
                exit $exitcode
            fi
            '''

    rule Barcode_15_CRISPResso2_ABE_list:
        input:
            expand("results/{prefix}_Barcode_pair_Results/ABE_Alleles_frequency_subgenome/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"]) if config["params"]["differentiation_subgenomes"] else expand("results/{prefix}_Barcode_pair_Results/ABE_Alleles_frequency/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"])
        output:
            "results/{prefix}_Barcode_pair_Results/All_ABE_samples_Barcode.txt"
        shell:
            '''
            cat {input} > {output}
            '''

if config["module"]["Barcode"]["CBE"]:
    rule Barcode_16_CRISPResso2_CBE_command:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            command="CRISPResso -r1 {params.workdir}/{input.fq} --output_folder {params.workdir}/{params.outdir}/CBE --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from C --conversion_nuc_to T {params.extr} "
            sed "s|^|${{command}} |" {input.description} > {output.bash_command}
            '''

    rule Barcode_17_CRISPResso2_CBE_command_subgenome:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE_subgenome.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            for i in `ls {input.fq} | sort -V | xargs`
            do
            subgenome=`echo ${{i}} | sed "s/\\.fastq\\.gz//"`
            command="CRISPResso -r1 {params.workdir}/{input.fq}/${{i}} --output_folder {params.workdir}/{params.outdir}/CBE_subgenome --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from C --conversion_nuc_to T {params.extr} "
            cat {input.description} | awk -v a="${{subgenome}}" -v command="${{command}}" '$2==a {{print command,$0}}' >> {output.bash_command}
            done
            '''

    rule Barcode_18_CRISPResso2_CBE:
        input:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE_subgenome.sh" if config["params"]["differentiation_subgenomes"] else "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_CBE.sh",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            CRISPResso_png = report(directory("results/{prefix}_Barcode_pair_Results/CBE_subgenome/CRISPResso_on_{Barcode_sample}_subgenomeA"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE figure") if config["params"]["differentiation_subgenomes"] else report(directory("results/{prefix}_Barcode_pair_Results/CBE/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE figure"),
            CRISPResso = report("results/{{prefix}}_Barcode_pair_Results/CBE_Alleles_frequency_subgenome/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE table") if config["params"]["differentiation_subgenomes"] else report("results/{{prefix}}_Barcode_pair_Results/CBE_Alleles_frequency/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_CBE.rst", category="2-1. Barcode CBE table")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            zaiti = config["zaiti"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_CRISPResso2_CBE_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.workdir}/{output.CRISPResso_png}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.CRISPResso}
                exit 0
            else
                echo "Success run {output.CRISPResso}" >> {params.workdir}/{log} 2>&1
                rows=`cat {input.description} | cut -f 6 -d " " | wc -l`
                for (( i=1; i<=$rows; i++ ))
                do
                    sgRNA=`cat {input.description} | cut -f 6 -d " " | sed -n ${{i}}p`
                    cat {params.workdir}/{output.CRISPResso_png}/Alleles_frequency_table_around_sgRNA_${{sgRNA}}.txt | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $3 ~ /False/ && $8 >= 0.2 {{print $0}}' | awk -F"\\t" -v a="{params.sampleid}" -v b="{params.prefix}" -v c="{params.zaiti}" -v d="sgRNA${{i}}" '{{SUM+=$8}}END{{print b"_"a"\\t"SUM"\\t"c"\\t"d}}' >> {output.CRISPResso}
                done
                exit $exitcode
            fi
            '''

    rule Barcode_19_CRISPResso2_CBE_list:
        input:
            expand("results/{prefix}_Barcode_pair_Results/CBE_Alleles_frequency_subgenome/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"]) if config["params"]["differentiation_subgenomes"] else expand("results/{prefix}_Barcode_pair_Results/CBE_Alleles_frequency/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"])
        output:
            "results/{prefix}_Barcode_pair_Results/All_CBE_samples_Barcode.txt"
        shell:
            '''
            cat {input} > {output}
            '''

if config["module"]["Barcode"]["cpf1"]:
    rule Barcode_20_CRISPResso2_cpf1_command:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            command="CRISPResso -r1 {params.workdir}/{input.fq} --output_folder {params.workdir}/{params.outdir}/cpf1 {params.extr} --cleavage_offset -6 --plot_window_size 25 --default_min_aln_score 50 "
            sed "s|^|${{command}} |" {input.description} > {output.bash_command}
            '''

    rule Barcode_21_CRISPResso2_cpf1_command_subgenome:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1_subgenome.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            for i in `ls {input.fq} | sort -V | xargs`
            do
            subgenome=`echo ${{i}} | sed "s/\\.fastq\\.gz//"`
            command="CRISPResso -r1 {params.workdir}/{input.fq}/${{i}} --output_folder {params.workdir}/{params.outdir}/cpf1_subgenome {params.extr} --cleavage_offset -6 --plot_window_size 25 --default_min_aln_score 50 "
            cat {input.description} | awk -v a="${{subgenome}}" -v command="${{command}}" '$2==a {{print command,$0}}' >> {output.bash_command}
            done
            '''

    rule Barcode_22_CRISPResso2_cpf1:
        input:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1_subgenome.sh" if config["params"]["differentiation_subgenomes"] else "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_cpf1.sh",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            CRISPResso_png = report(directory("results/{prefix}_Barcode_pair_Results/cpf1_subgenome/CRISPResso_on_{Barcode_sample}_subgenomeA"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_cpf1.rst", category="2-1. Barcode cpf1 figure") if config["params"]["differentiation_subgenomes"] else report(directory("results/{prefix}_Barcode_pair_Results/cpf1/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_cpf1.rst", category="2-1. Barcode cpf1 figure"),
            CRISPResso = report("results/{{prefix}}_Barcode_pair_Results/cpf1_Alleles_frequency_subgenome/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_cpf1.rst", category="2-1. Barcode cpf1 table") if config["params"]["differentiation_subgenomes"] else report("results/{{prefix}}_Barcode_pair_Results/cpf1_Alleles_frequency/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_cpf1.rst", category="2-2. Barcode cpf1 table")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            zaiti = config["zaiti"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_CRISPResso2_cpf1_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.workdir}/{output.CRISPResso_png}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.CRISPResso}
                exit 0
            else
                echo "Success run {output.CRISPResso}" >> {params.workdir}/{log} 2>&1
                rows=`cat {input.description} | cut -f 6 -d " " | wc -l`
                for (( i=1; i<=$rows; i++ ))
                do
                    sgRNA=`cat {input.description} | cut -f 6 -d " " | sed -n ${{i}}p`
                    cat {params.workdir}/{output.CRISPResso_png}/Alleles_frequency_table_around_sgRNA_${{sgRNA}}.txt | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $3 ~ /False/ && $8 >= 0.2 {{print $0}}' | awk -F"\\t" -v a="{params.sampleid}" -v b="{params.prefix}" -v c="{params.zaiti}" -v d="sgRNA${{i}}" '{{SUM+=$8}}END{{print b"_"a"\\t"SUM"\\t"c"\\t"d}}' >> {output.CRISPResso}
                done
                exit $exitcode
            fi
            '''

    rule Barcode_23_CRISPResso2_cpf1_list:
        input:
            expand("results/{prefix}_Barcode_pair_Results/cpf1_Alleles_frequency_subgenome/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"]) if config["params"]["differentiation_subgenomes"] else expand("results/{prefix}_Barcode_pair_Results/cpf1_Alleles_frequency/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], Barcode_sample = Barcode_samples, zaiti = config["zaiti"])
        output:
            "results/{prefix}_Barcode_pair_Results/All_cpf1_samples_Barcode.txt"
        shell:
            '''
            cat {input} > {output}
            '''

if config["module"]["Barcode"]["PE"]:
    rule Barcode_24_CRISPResso2_PE_command:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}.fastq.gz",
            description = "results/{prefix}_Barcode_pair_Results/description_PE/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_PE.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            command="CRISPResso -r1 {params.workdir}/{input.fq} --output_folder {params.workdir}/{params.outdir}/PE {params.extr} "
            sed "s|^|${{command}} |" {input.description} > {output.bash_command}
            '''

    rule Barcode_25_CRISPResso2_PE_command_subgenome:
        input:
            fq = "results/{prefix}_Barcode_pair_Results/Split_fastqs/{Barcode_sample}_subgenome_fastq",
            description = "results/{prefix}_Barcode_pair_Results/description/{Barcode_sample}_description.txt"
        output:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_PE_subgenome.sh"
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            outdir = lambda w, output: "/".join(os.path.dirname(output.bash_command).split("/")[0:2]),
            extr = config["params"]["CRISPResso"]
        shell:
            '''
            for i in `ls {input.fq} | sort -V | xargs`
            do
            subgenome=`echo ${{i}} | sed "s/\\.fastq\\.gz//"`
            command="CRISPResso -r1 {params.workdir}/{input.fq}/${{i}} --output_folder {params.workdir}/{params.outdir}/PE_subgenome {params.extr} "
            cat {input.description} | awk -v a="${{subgenome}}" -v command="${{command}}" '$2==a {{print command,$0}}' >> {output.bash_command}
            done
            '''

    rule Barcode_26_CRISPResso2_PE:
        input:
            bash_command = "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_PE_subgenome.sh" if config["params"]["differentiation_subgenomes"] else "results/{prefix}_Barcode_pair_Results/Bash_Command/{Barcode_sample}_CRISPResso2_PE.sh",
            description = "results/{prefix}_Barcode_pair_Results/description_PE/{Barcode_sample}_description.txt"
        output:
            CRISPResso_png = report(directory("results/{prefix}_Barcode_pair_Results/PE_subgenome/CRISPResso_on_{Barcode_sample}_subgenomeA"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_PE.rst", category="2-1. Barcode PE figure") if config["params"]["differentiation_subgenomes"] else report(directory("results/{prefix}_Barcode_pair_Results/PE/CRISPResso_on_{Barcode_sample}"), patterns=["9.Alleles_frequency_table_around_sgRNA_{name}.png"], caption="../report/Barcode_PE.rst", category="2-1. Barcode PE figure"),
            CRISPResso = report("results/{{prefix}}_Barcode_pair_Results/PE_Alleles_frequency_subgenome/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_PE.rst", category="2-1. Barcode PE table") if config["params"]["differentiation_subgenomes"] else report("results/{{prefix}}_Barcode_pair_Results/PE_Alleles_frequency/{{Barcode_sample}}_{{prefix}}_{}_All_Alleles_frequency.txt".format(config["zaiti"]), caption="../report/Barcode_PE.rst", category="2-1. Barcode PE table")
        threads: 1
        params:
            workdir = config["workdir"],
            prefix = config["prefix"],
            sampleid = "{Barcode_sample}",
            zaiti = config["zaiti"]
        conda:
            "../envs/Barcode_Hi-Tom.yaml"
        log:
            "logs/{prefix}/Barcode_CRISPResso2_PE_{prefix}-{Barcode_sample}.log"
        shell:
            '''
            set +e
            bash {input.bash_command} > {log} 2>&1
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
                touch {params.workdir}/{output.CRISPResso_png}/9.Alleles_frequency_table_around_sgRNA_ATCG.png >> {params.workdir}/{log} 2>&1
                touch {output.CRISPResso}
                exit 0
            else
                echo "Success run {output.CRISPResso}" >> {params.workdir}/{log} 2>&1
                rows=`cat {input.description} | cut -f 6 -d " " | wc -l`
                for (( i=1; i<=$rows; i++ ))
                do
                    sgRNA=`cat {input.description} | cut -f 6 -d " " | sed -n ${{i}}p`
                    cat {params.workdir}/{output.CRISPResso_png}/Prime-edited.Alleles_frequency_table_around_PE-Extension.txt | awk -F"\\t" 'BEGIN{{OFS="\\t"}} $3 ~ /False/ && $8 >= 0.2 {{print $0}}' | awk -F"\\t" -v a="{params.sampleid}" -v b="{params.prefix}" -v c="{params.zaiti}" -v d="sgRNA${{i}}" '{{SUM+=$8}}END{{print b"_"a"\\t"SUM"\\t"c"\\t"d}}' >> {output.CRISPResso}
                done
                exit $exitcode
            fi
            '''

    rule Barcode_27_CRISPResso2_PE_list:
        input:
            expand("results/{prefix}_Barcode_pair_Results/PE_Alleles_frequency_subgenome/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], zaiti = config["zaiti"], Barcode_sample = Barcode_samples) if config["params"]["differentiation_subgenomes"] else expand("results/{prefix}_Barcode_pair_Results/PE_Alleles_frequency/{Barcode_sample}_{prefix}_{zaiti}_All_Alleles_frequency.txt", prefix = config["prefix"], zaiti = config["zaiti"], Barcode_sample = Barcode_samples)
        output:
            "results/{}_Barcode_pair_Results/All_PE_samples_Barcode.txt".format(config["prefix"])
        shell:
            '''
            cat {input} > {output}
            '''
