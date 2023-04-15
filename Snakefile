import yaml
from pathlib import Path

BASE_DIR = Path(workflow.basedir)

### Validation of schemas ###
##### load config and sample sheets ##
configfile: BASE_DIR / "config/config.yaml"

##### load rules #####
include: "rules/common.smk"  # python helper functions

paths = create_path_accessor()

hifi_path = expand(paths.hifiasm.hifiasm + "/" + config["prefix"] + ".asm.dip.{hap}.p_ctg.gfa", hap = ["hap1", "hap2"])

# rule test:
#     input:
#         config["hifi"]
#     shell:
#         "echo {paths}"


rule all:
    input:
        # "results/wtdbg2/father/father.cns.fa"
        expand("results/merge/{hap}/output/scaffold.fasta", hap=["mother", "father"])

rule yak_father:
    input:
        config["father_illumina"],
    output:
        paths.yak.father
    threads: 16
    shell:
        "yak count -k31 -b37 -t {threads} -o {output} {input}"

rule yak_mother:
    input:
        config["mother_illumina"],
    output:
        paths.yak.mother
    threads: 16
    shell:
        "yak count -k31 -b37 -t {threads} -o {output} {input}"

rule hifiasm:
    input:
        hifi = config["hifi"],
        patyak = paths.yak.mother,
        matyak = paths.yak.father
    output: 
        hap1 = paths.hifiasm.hap1,
        hap2 = paths.hifiasm.hap2,
    params: prefix = config["prefix"],
            path = config["output_dir"]
    threads: 32
    shell:
        '''
        hifiasm -o {params.prefix} -t {threads} -1 {input.patyak} -2 {input.matyak} {input.hifi}
        mkdir -p {params.path}/hifiasm
        mv {params.prefix}* {params.path}/hifiasm
        '''

rule gfa_tofa:
    input: 
        hap1 = paths.hifiasm.hap1,
        hap2 = paths.hifiasm.hap2,
    output: 
        hap1fa = paths.hifiasm.hap1fa,
        hap2fa = paths.hifiasm.hap2fa,
    shell:
        '''
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1} > {output.hap1fa}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2} > {output.hap2fa}
        '''

rule haplotype:
    input:
        reads = config["nanopore"],
        father = config["father_illumina"],
        mother = config["mother_illumina"]
    output: 
        Fat = "results/canu/haplotype/haplotype-Fat.fasta.gz",
        Mat = "results/canu/haplotype/haplotype-Mat.fasta.gz",
        unknown = "results/haplotype/haplotype-unknown.fasta.gz"
    params: 
        path = "results/canu/",
        outbase = config['basename'] + ".hifiasm_trio",
        readsType = config["readsType"],
        prefix = config["prefix"],
        genomeSize = config["genomeSize"],
    threads: 32
    shell:
        '''
        mkdir -p {params.path}
        canu -haplotype -p {params.path}{params.prefix} -d haplotype genomeSize={params.genomeSize} 
        -haplotypeFather {input.father} 
        -haplotypeMother {input.mother} 
        -{params.readsType} {input.reads}
        rm -rf {params.path}haplotype/0-kmers
        '''

rule wtdbg2:
    input:
        Fat = "results/canu/haplotype/haplotype-Fat.fasta.gz",
        Mat = "results/canu/haplotype/haplotype-Mat.fasta.gz",
        unknown = "results/haplotype/haplotype-unknown.fasta.gz"
    params: 
        path = "results/wtdbg2/",
        readsType = "ont",
        genomeSize = "3.1g"
    output: 
            "results/wtdbg2/father/father.cns.fa",
            "results/wtdbg2/mother/mother.cns.fa",
            touch("wtdbg2.done")
    threads: 32
    shell:
        """
        wtdbg2.pl -t {threads} -x {params.readsType} -g {params.genomeSize} -o {params.path}father/father {input.Fat} {input.unknown}

        wtdbg2.pl -t {threads} -x {params.readsType} -g {params.genomeSize} -o {params.path}mother/mother {input.Mat} {input.unknown}

        echo {BASE_DIR}/{input.Fat} > lgs_father.fofn
        echo {BASE_DIR}/{input.Mat} > lgs_mother.fofn
        echo {BASE_DIR}/{input.unknown} >> *.fofn
        """

rule polish:
    input:
        "wtdbg2.done"
    params: 
        father_cfg = config["father_cfg"],
        mother_cfg = config["mother_cfg"],
        path = "results/polish/",
    output: 
            "results/polish/father/genome.nextpolish.fasta",
            "results/polish/mother/genome.nextpolish.fasta",
            touch("polish.done")
    threads: 32
    shell:
        """
        BD={params.path}
        mkdir -p $BD"father" $BD"mother"
        cp {params.father_cfg} $BD
        cp {params.mother_cfg} $BD
        mv *.fofn $BD
        cd $BD

        nextDenovo father.cfg
        nextDenovo mother.cfg

        echo {BASE_DIR}
        """

rule merge:
    input:
        # "polish.done",
        hifi_father = paths.hifiasm.hap1fa,
        hifi_mother = paths.hifiasm.hap2fa,
        ont_father = "results/polish/father/genome.nextpolish.fasta",
        ont_mother = "results/polish/mother/genome.nextpolish.fasta",
    output: 
        expand("results/merge/{hap}/output/scaffold.fasta", hap=["mother", "father"])
    params: 
        prefix = "NA12878",
        genomeSize = "3.1g",
    threads: 32
    shell:
        """
        mkdir -p results/merge/{{father,mother}}/{{input,output,temp}}
        cd results/merge
        wget https://github.com/bioinfomaticsCSU/MAC/archive/refs/tags/V2.1.zip
        unzip *.zip
        mv MAC*/*cpp ./
        mv *cpp mac.cpp
        rm -rf MAC* *.zip
        g++ mac.cpp -o MAC

        cp {input.hifi_father} father/input/hap1.fa
        cp {input.ont_father} father/input/father.fa
        cp {input.hifi_mother} mother/input/hap2.fa
        cp {input.ont_mother} mother/input/mother.fa
        
        cd father
        ../mac hap1.fa father.fa

        cd ../mother
        ../mac hap2.fa mother.fa

        """
        # make sure the path of MUMmer has been added into system variable.


rule scafflod:
    input:
        reads = "",
        father = "",
        mother = ""
    output: "assembly/nanopore/"
    params: 
        readsType = "nanopore",
        prefix = "NA12878",
        genomeSize = "3.1g",
    threads: 32
    shell:
        """
        echo 1
        """