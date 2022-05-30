import os
import re
import yaml
import pathlib
from typing import List
import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.15")


#======================================================
###### cluster directory
#======================================================
logs = pathlib.Path("./logs/cluster")

if logs.is_dir():
    pass
    #print("#"*40 + "\n" + "        The logs exist." + "\n" + "#"*40 + "\n")
else:
    #print("#"*40 + "\n" + "Now, create the logs directory." + "\n")
    os.system('mkdir -p ./logs/cluster')
    #print("#"*40 + "\n")

#======================================================
###### Config files and sample sheets
#======================================================

configfile: "config.yaml"
validate(config, schema="./schemas/config.schema.yaml")

Pfam_ids = pd.read_table(config["domain_metadata"]).set_index("Domain_ID", drop=False)
validate(Pfam_ids, schema="schemas/domain_metadata.schema.yaml")

#======================================================
##### Wildcard constraints
#======================================================

#======================================================
##### Helper functions
#======================================================

#======================================================
###### target rules
#======================================================

rule all:
    input:
        "Gene_Family/All_{}_pep.fa".format(config["prefix"]),
        "Gene_Family/All_{}_CC_pep.fa".format(config["prefix"]) if config["module"]["paircoil2"] else []
            


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"

## select which rule to run
include: "rules/Identify_Gene_Family.smk"
