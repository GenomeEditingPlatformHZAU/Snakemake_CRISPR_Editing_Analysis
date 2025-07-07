{% if snakemake.config["params"]["CRISPResso"] %}
ABE单碱基分析结果，其中额外参数 ``{{ snakemake.config["params"]["CRISPResso"] }}`` 运行CRISPResso.
{% else %}
ABE单碱基分析结果，其中默认参数运行CRISPResso。
{% endif %}

仅显示大于0.2%的编辑统计结果。
