# Snakemake workflow: CRISPR sgRNA Designer

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html)

This workflow performs gene family identification.

## Authors

* Zhongping Xu (@hopetogy), http://tiramisutes.github.io/

## Usage

### Simple

#### Step 1: Install workflow / git clone

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

Git clone this repository using follows command 👇

    git clone --recursive https://github.com/tiramisutes/Snakemake_Gene_Family.git

目录结构如下：

```
|-- bin             # 相关软件目录
|-- cluster.json    # 集群运行相关设置
|-- config.yaml     # 本snakemake控制文件
|-- envs            # rule运行环境设置
|-- LICENSE         # 协议
|-- README.md       # 说明
|-- report          # html报告注释
|-- Resources       # 运行前需准备资源文件夹
|-- rules           # 运行rule文件夹
|-- schemas         # config.yaml 文件源文件
|-- scripts         # 本snakemake 用到相关额外脚本
`-- Snakefile       # 控制rule运行顺序等
```

#### Step 2: Prepare data and software

#### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

or

    snakemake --use-conda --jobs 100 --cluster-config cluster.json --cluster "bsub -q {cluster.queue} -o {cluster.output} -e {cluster.error}"

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

#### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html).
