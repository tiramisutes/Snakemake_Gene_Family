rule InterProScan_install:
    output:
        "bin/{}/InterProScan_install.log".format(config["InterProScan_name"])
    threads: 1
    conda:
        "../envs/Identify_Gene_Family.yaml"
    params:
        workdir = config["workdir"],
        InterProScan_url = config["InterProScan_url"],
        InterProScan_name = config["InterProScan_name"]
    shell:
        '''
        cd {params.workdir}/bin
        inproname=`basename {params.InterProScan_url}`
        #wget {params.InterProScan_url}
        #tar -pxvzf ${{inproname}}
        cd {params.InterProScan_name}
        python3 initial_setup.py > {params.workdir}/{output} 2>&1
        '''

def genome_pep_name():
    import os
    genome_pep = config["genome_pep"]
    return os.path.basename(genome_pep)

if config["module"]["InterProScan"]:
    rule InterProScan_genome:
        input:
            genome_pep = config["genome_pep"] if config["genome_pep"] else "Resources/{}_pep.fa".format(config["prefix"]),
            inst = "bin/{}/InterProScan_install.log".format(config["InterProScan_name"])
        output:
            "InterProScan/{}.tsv".format(genome_pep_name())
        threads: 40
        params:
            workdir = config["workdir"],
            InterProScan_name = config["InterProScan_name"],
            outdir = lambda w, output: os.path.dirname(output[0])
        conda:
            "../envs/Identify_Gene_Family.yaml"
        log:
            "logs/InterProScan_genome.log"
        shell:
            '''
            cd {params.workdir}/bin/{params.InterProScan_name}
            ./interproscan.sh -appl Pfam -i {input.genome_pep} \
            -goterms -iprlookup -pa -f TSV,XML,GFF3,HTML \
            -d {params.workdir}/{params.outdir} \
            -cpu {threads} > {params.workdir}/{log} 2>&1
            '''

rule Download_Pfam:
    output:
        "Resources/Pfam_hmm/{Domain_ID}.hmm"
    threads: 1
    params:
        workdir = config["workdir"],
        Pfam_id = "{Domain_ID}",
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/Download_Pfam_{Domain_ID}.log"
    shell:
        '''
        mkdir -p {params.outdir}; cd {params.outdir}
        wget http://pfam.xfam.org/family/{params.Pfam_id}/hmm -O {params.Pfam_id}.hmm > {params.workdir}/{log} 2>&1

        '''

rule hmmsearch_first:
    input:
        genome_pep = config["genome_pep"] if config["genome_pep"] else "Resources/{}_pep.fa".format(config["prefix"]),
        hmm = "Resources/Pfam_hmm/{Domain_ID}.hmm"
    output:
        "HMMER/{Domain_ID}_Results.txt"
    threads: 20
    conda:
        "../envs/Identify_Gene_Family.yaml"
    log:
        "logs/hmmsearch_first_{Domain_ID}.log"
    shell:
        '''
        hmmsearch --cut_tc --cpu {threads} --domtblout {output} {input.hmm} {input.genome_pep} > {log} 2>&1
        '''

rule hmmsearch_first_fa:
    input:
        genome_pep = config["genome_pep"] if config["genome_pep"] else "Resources/{}_pep.fa".format(config["prefix"]),
        hmm1 = "HMMER/{Domain_ID}_Results.txt"
    output:
        qua_id = "HMMER/{Domain_ID}_qua_id.txt",
        qua_fa = "HMMER/{Domain_ID}_qua.fa"
    threads: 1
    conda:
        "../envs/Identify_Gene_Family.yaml"
    params:
        evalue = config["hmmsearch"]["fir_evalue"],
        outdir = lambda w, output: os.path.dirname(output.qua_id)
    log:
        "logs/hmmsearch_first_fa_{Domain_ID}.log"
    shell:
        '''
        cat {input.hmm1} | grep -v "#" | awk -v e="{params.evalue}" '($7 + 0) < e' | cut -f1 -d  " " | sort -u > {output.qua_id}
        if [[ ! -s {output.qua_id} ]]
        then
        echo "This {output.qua_id} is empty when filter use {params.evalue} e-value for hmmsearch (first run) results. So, no filtes." >> {params.outdir}/note.txt
        cat {input.hmm1} | grep -v "#" | cut -f1 -d  " " | sort -u > {output.qua_id}
        fi
        seqtk subseq {input.genome_pep} {output.qua_id} > {output.qua_fa}
        '''

rule hmmsearch_first_fa_clustalw:
    input:
        qua_fa = "HMMER/{Domain_ID}_qua.fa"
    output:
        qua_aln = "HMMER/{Domain_ID}_qua.aln",
        qua_dnd = "HMMER/{Domain_ID}_qua.dnd"
    threads: 1
    conda:
        "../envs/Identify_Gene_Family.yaml"
    log:
        "logs/hmmsearch_first_fa_clustalw_{Domain_ID}.log"
    shell:
        '''
        echo -e "1\\n{input.qua_fa}\\n2\\n1\\n{output.qua_aln}\\n{output.qua_dnd}\\nX\\n\\nX\\n" | clustalw2 > {log} 2>&1
        '''

rule hmmsearch_second:
    input:
        genome_pep = config["genome_pep"] if config["genome_pep"] else "Resources/{}_pep.fa".format(config["prefix"]),
        qua_aln = "HMMER/{Domain_ID}_qua.aln"
    output:
        qua_hmm = "HMMER/{Domain_ID}_qua.hmm",
        qua_result = "HMMER/{Domain_ID}_qua_Results.txt",
        family_hmm = "HMMER/{Domain_ID}_hmm_list.txt"
    threads: 20
    conda:
        "../envs/Identify_Gene_Family.yaml"
    params:
        evalue = config["hmmsearch"]["se_evalue"],
        outdir = lambda w, output: os.path.dirname(output.qua_hmm)
    log:
        "logs/hmmsearch_second_{Domain_ID}.log"
    shell:
        '''
        hmmbuild {output.qua_hmm} {input.qua_aln} >> {log} 2>&1
        hmmsearch --cpu {threads} --domtblout {output.qua_result} {output.qua_hmm} {input.genome_pep} >> {log} 2>&1
        cat {output.qua_result} | grep -v "#" | awk -v e="{params.evalue}" '($7 + 0) < e' > {output.family_hmm}
        if [[ ! -s {output.family_hmm} ]]
        then
        echo "This {output.family_hmm} is empty when filter use {params.evalue} e-value for hmmsearch (second run) results. So, no filtes." >> {params.outdir}/note.txt
        cat {output.qua_result} | grep -v "#" | > {output.family_hmm}
        fi
        '''

if config["module"]["InterProScan"]:
    pass
else:
    rule InterProScan_run:
        input:
            gene_pep = "HMMER/{Domain_ID}_pep.fa",
            inst = "bin/{}/InterProScan_install.log".format(config["InterProScan_name"])
        output:
            "InterProScan/{Domain_ID}_pep.fa.tsv"
        threads: 40
        params:
            workdir = config["workdir"],
            InterProScan_name = config["InterProScan_name"],
            outdir = lambda w, output: os.path.dirname(output[0])
        conda:
            "../envs/Identify_Gene_Family.yaml"
        log:
            "logs/InterProScan_run_{Domain_ID}.log"
        shell:
            '''
            cd {params.workdir}/bin/{params.InterProScan_name}
            ./interproscan.sh -appl Pfam -i {params.workdir}/{input.gene_pep} \
                -goterms -iprlookup -pa -f TSV,XML,GFF3,HTML \
                -d {params.workdir}/{params.outdir} -cpu {threads} > {params.workdir}/{log} 2>&1
            '''

    rule InterProScan_run_merge:
        input:
            expand("InterProScan/{u.Domain_ID}_pep.fa.tsv", u = Pfam_ids.itertuples())
        output:
            report("Gene_Family/All_{}_pep.fa.tsv".format(config["prefix"]), caption="../report/InterProScan_run_merge.rst", category="3. InterProScan annotation of all gene family")
        threads: 1
        log:
            "logs/InterProScan_run_merge.log"
        shell:
            '''
            cat {input} > {output}
            '''

rule family_fasta:
    input:
        genome_pep = config["genome_pep"] if config["genome_pep"] else "Resources/{}_pep.fa".format(config["prefix"]),
        family_hmm = "HMMER/{Domain_ID}_hmm_list.txt",
        family_interproscan = "InterProScan/{}.tsv".format(genome_pep_name()) if config["module"]["InterProScan"] else "Gene_Family/All_{}_pep.fa.tsv".format(config["prefix"])
    output:
        hmm_id = "HMMER/{Domain_ID}_hmm_id.txt",
        interproscan_id = "HMMER/{Domain_ID}_interproscan_id.txt" if config["module"]["InterProScan"] else [],
        gene_id = "HMMER/{Domain_ID}_id.txt",
        gene_pep = report("HMMER/{Domain_ID}_pep.fa", caption="../report/Gene_family_pep.rst", category="1. Protein sequences of gene family for each Pfam ID.")
    threads: 1
    conda:
        "../envs/Identify_Gene_Family.yaml"
    params:
        Pfam_id = "{Domain_ID}"
    log:
        "logs/family_fasta_{Domain_ID}.log"
    shell:
        '''
        cat {input.family_hmm} | cut -f1 -d  " " | sort -u > {output.hmm_id}
        family_interproscan="{input.family_interproscan}"
        interproscan_id="{output.interproscan_id}"
        if [ ! -n "${{family_interproscan}}" ]; then
            echo "No interproscan input."
        else
            cat ${{family_interproscan}} | \\grep -w {params.Pfam_id} | cut -f1 | sort -u > ${{interproscan_id}}
        fi
        cat {output.hmm_id} {output.interproscan_id} | sort -u > {output.gene_id}
        seqtk subseq {input.genome_pep} {output.gene_id} > {output.gene_pep}
        '''

if config["module"]["paircoil2"]:
    rule paircoil2:
        input:
            "HMMER/{Domain_ID}_pep.fa"
        output:
            "paircoil2/{Domain_ID}_CC_pep.fa"
        threads: 1
        conda:
            "../envs/Identify_Gene_Family.yaml"
        params:
            workdir = config["workdir"],
            Pfam_id = "{Domain_ID}"
        log:
            "logs/paircoil2_{Domain_ID}.log"
        shell:
            '''
            cd {params.workdir}/bin/paircoil2
            ./paircoil2 -win 28 {params.workdir}/{input} {params.workdir}/{output} > {params.workdir}/{log} 2>&1
            '''

    rule merge_family_fasta:
        input:
            fa = expand("HMMER/{u.Domain_ID}_pep.fa", u = Pfam_ids.itertuples()),
            cc_fa = expand("paircoil2/{u.Domain_ID}_CC_pep.fa", u = Pfam_ids.itertuples())
        output:
            fa = report("Gene_Family/All_{}_pep.fa".format(config["prefix"]), caption="../report/merge_family_fasta.rst", category="2. Protein sequences of all gene family"),
            cc_fa = "Gene_Family/All_{}_CC_pep.fa".format(config["prefix"])
        threads: 1
        log:
            "logs/merge_family_fasta.log"
        shell:
            '''
            cat {input.fa} > {output.fa}
            cat {input.cc_fa} > {output.cc_fa}
            '''
else:
    rule merge_family_fasta:
        input:
            fa = expand("HMMER/{u.Domain_ID}_pep.fa", u = Pfam_ids.itertuples())
        output:
            fa = report("Gene_Family/All_{}_pep.fa".format(config["prefix"]), caption="../report/merge_family_fasta.rst", category="2. Protein sequences of all gene family")
        threads: 1
        log:
            "logs/merge_family_fasta.log"
        shell:
            '''
            cat {input.fa} > {output.fa}
            '''
