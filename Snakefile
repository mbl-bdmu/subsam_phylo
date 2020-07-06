configfile: "config.yaml"

rule all:
    input:
        expand(
        	"{folders}/{sample}.{type}", 
        	folders=config["folders"], 
        	sample=config["sample"], 
        	type=config["types"]
        	)

# Rules for generating plots & output files.
rule extract_sequences:
    input:
        "master/{sample}.fa",
        "{folder}/{sample}.headers.txt"
    output:
        "{folder}/{sample}.fa"
    script:
        "scripts/io_alignment.py"

rule generate_fasttree:
    input:
        "{folder}/{sample}.fa"
    output:
        "{folder}/{sample}.fasttre.nwk"
    shell:
        "scripts/run_fasttree.sh {input}"

rule root_to_tip:
    input:
        "{folder}/{sample}.headers.txt",
        "{folder}/{sample}.fasttre.nwk",
        "{folder}/{sample}.fa"
    output:
        "{folder}/{sample}.dates.txt",
        "{folder}/{sample}.rtt.png" 
    shell:
        "scripts/plot_rtt.sh {input} {output}"

rule tree_img:
    input:
        "{folder}/{sample}.fasttre.nwk"
    output:
        "{folder}/{sample}.fasttre.nwk.png",
        "{folder}/{sample}.fasttre.nwk.pdf"
    shell:
        "Rscript --vanilla scripts/plot_tree.R {input}"

# Rules for subsampling schemes.
rule sub_random:
    input:
        "sub/{sample}.headers.txt"
    output:
        "sub_random/{sample}.headers.txt"
    params:
        size=config["randsize"]
    script:
        "scripts/subsam_random.py"

rule sub_rtl:
    input:
        "sub/{sample}.fasttre.nwk"
    output:
        "sub_rtl/{sample}.headers.txt"
    params:
        rtl=config["rtl"],
        threads=config["rtl_threads"]
    shell:
        "scripts/subsam_rtl.sh {input} {output} {params.rtl} {params.threads}"

rule sub_uniform_tempo:
    input:
        "sub/{sample}.headers.txt"
    output:
        "sub_uniform-tempo/{sample}.headers.txt"
    params:
        size=config["unisize"],
        intervals=config["uniinterv"]
    script:
        "scripts/subsam_uniform_tempo.py"

rule sub_uniform_spatiotempo:
    input:
        "sub/{sample}.headers.txt",
        "master/{sample}.loc.tsv"
    output:
        "sub_uniform-spatiotempo/{sample}.headers.txt"
    params:
        size=config["unisize"],
        intervals=config["uniinterv"]
    script:
        "scripts/subsam_uniform_spatiotempo.py"

rule sub_uniform_spatiotempo_alt:
    input:
      "sub/{sample}.headers.txt",
      "master/{sample}.loc.tsv"
    output:
      "sub_uniform-spatiotempo-alt/{sample}.headers.txt"
    script:
      "scripts/subsam_uniform_spatiotempo_alt.py"

rule sub_postsubsam:
    input:
        "sub/{sample}.headers.txt",
        "master/{sample}.loc.tsv"
    output:
        "sub_postsubsampling/{sample}.headers.txt"
    params:
        size=config["fixperc"]
    script:
        "scripts/subsam_post.py"
