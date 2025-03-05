rule Barcode_Design:
    output:
        barcode = report("results/Barcode_Primer/{}_barcode.txt".format(config["prefix"]), caption="../report/Barcode_Design.rst", category="1-1. Barcode Design Raw"),
        ERbarcode = report("results/Barcode_Primer/ER_{}_barcode.txt".format(config["prefix"]), caption="../report/Barcode_Design_company.rst", category="1-2. Barcode Design Company"),
        ERObarcode = report("results/Barcode_Primer/ERO_{}_barcode.txt".format(config["prefix"]), caption="../report/Barcode_Design_Oneline.rst", category="1-3. Barcode Design Oneline"),
        ERTMbarcode = report("results/Barcode_Primer/ER_TM_{}_barcode.txt".format(config["prefix"]), caption="../report/Barcode_Design_TM.rst", category="1-4. Barcode Design Tm")
    threads: 1
    conda:
        "../envs/Barcode_Hi-Tom.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        primerF = config["params"]["primerF"],
        primerR = config["params"]["primerR"],
        length = config["params"]["length"],
        total = config["params"]["total"],
        minimum_dif = config["params"]["minimum_dif"]
    log:
        "logs/{}/Barcode_Design.log".format(config["prefix"])
    shell:
        '''
        python {params.workdir}/workflow/bin/BarCode/barcode_generator.py \
         --project {params.prefix} \
         --primerF {params.primerF} \
         --primerR {params.primerR} \
         --length {params.length} \
         --total {params.total} \
         --minimum_dif {params.minimum_dif} \
         --output {output.barcode} >> {log} 2>&1
        
        cat {output.barcode} | \\grep "{params.prefix}" | sed 'N;s/\\n/\\t/' | cat -n | awk 'BEGIN{{print "Primer_name,Seq";OFS="\\t"}} {{print $1"_"$2,$3,$1"_"$4,$5}}' | sed "s/\\t/\\n/2" | sed "s/,/\\t/" > {output.ERbarcode}
        cat {output.barcode} | \\grep "{params.prefix}" | sed 'N;s/\\n/\\t/' | cat -n | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$4,$5}}' > {output.ERObarcode}
        perl {params.workdir}/workflow/bin/BarCode/getNearestNeighborTm.pl -c 2 -f {output.ERbarcode} > {output.ERTMbarcode}
        '''
