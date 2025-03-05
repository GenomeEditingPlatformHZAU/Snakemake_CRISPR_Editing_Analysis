{% if snakemake.config["params"]["CRISPResso"] %}
普通CRISPR/Cas9分析结果，其中CRISPResso运行额外参数为：``{{ snakemake.config["params"]["CRISPResso"] }}``。
{% else %}
普通CRISPR/Cas9分析结果，其中CRISPResso运行默认参数。
{% endif %}

仅显示大于0.2%的编辑统计结果。