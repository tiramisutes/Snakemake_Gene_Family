{% if snakemake.config["module"]["InterProScan"] %}
用 HMMER_ 软件经隐马可夫模型和 InterProScan_ 得到的基因家族成员。
{% else %}
用 HMMER_ 软件经隐马可夫模型得到的基因家族成员。
{% endif %}

.. _HMMER: http://hmmer.org/
.. _InterProScan: https://www.ebi.ac.uk/interpro/search/sequence/
