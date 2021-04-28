## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.yml"

rule all:
    input:
        "results/folder/result.ext"
        
# Rules

include: "rules/rule.smk"
