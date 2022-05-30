{% if snakemake.config["module"]["InterProScan"] %}
所有Pfam ID 用 HMMER_ 软件经隐马可夫模型和 InterProScan_ 进行基因家族成员的鉴定，并将所有结果进行合并。
{% else %}
所有Pfam ID 用 HMMER_ 软件经隐马可夫模型进行基因家族成员的鉴定，并将所有结果进行合并。
{% endif %}

.. _HMMER: http://hmmer.org/
.. _InterProScan: https://www.ebi.ac.uk/interpro/search/sequence/
