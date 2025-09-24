datasets_full = ['sts_drosophilia', 'MOB']
#, 'MF_AlzheimerBrain', 'MF_Hypothalamus'
for dataset in datasets_full:
    dataset_name = dataset
    datasets = [dataset_name, dataset_name + '/shuffle1', dataset_name + '/shuffle2', dataset_name + '/shuffle3', dataset_name + '/shuffle4', dataset_name + '/shuffle5']
alpha_in = 0.05


methods = [
    #'spagcn',
    #'scbsp',
    'smash',
    #'spagft',
    'spatialde',
    'spatialde2',
    #'spark',
    #'sparkx',
    #'spagene',
    'meringue'
]

rule all:
    input:
        expand("results/{dataset}/svg_counts.svg", dataset=datasets),
        #expand("results/{dataset}/spagcn_svgs.csv", dataset=datasets),
        #expand("results/{dataset}/scbsp_svgs.csv", dataset=datasets),
        expand("results/{dataset}/smash_svgs.csv", dataset=datasets),
        #expand("results/{dataset}/spagft_svgs.csv", dataset=datasets),
        expand("results/{dataset}/spatialde_svgs.csv", dataset=datasets),
        expand("results/{dataset}/spatialde2_svgs.csv", dataset=datasets),
        #expand("results/{dataset}/spark_svgs.csv", dataset=datasets),
        #expand("results/{dataset}/sparkx_svgs.csv", dataset=datasets),
        #expand("results/{dataset}/spagene_svgs.csv", dataset=datasets),
        expand("results/{dataset}/meringue_svgs.csv", dataset=datasets),
        expand("results/{dataset}/similarity_heatmap.svg", dataset=datasets)

rule comparison:
    input:
        method_outputs=lambda wildcards: [
            f"results/{wildcards.dataset}/{method}_svgs.csv" for method in methods
        ]
    output:
        svg_counts="results/{dataset}/svg_counts.svg",
        heatmap="results/{dataset}/similarity_heatmap.svg"
    params:
        methods=" ".join(methods),
        folder=lambda wildcards: f"results/{wildcards.dataset}",
        dataset=lambda wildcards: wildcards.dataset
    conda:
        "envs/spagcn_win.yaml"
    benchmark:
        "benchmarks/{dataset}/comparison.txt"
    shell:
        """
        python scripts/comparison.py \
            --inputs {params.methods} \
            --dataset {params.dataset} \
            --output {params.folder} \
            --folder {params.folder}
        """


rule spagcn:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spagcn_svgs.csv'
    params:
        n_neighbors = 30,
        alpha = alpha_in,
        dataset = '{dataset}'
    conda:
        'envs/spagcn_win.yaml'
    benchmark:
        "benchmarks/{dataset}/spagcn.txt"
    shell:
        """
        python scripts/spagcn.py \
            --alpha {params.alpha} \
            --coords {input.coords} \
            --counts {input.counts}\
            --dataset {params.dataset} \
            --n_neighbors 10 \
        """

rule scbsp:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/scbsp_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/r-env.yaml'
    benchmark:
        "benchmarks/{dataset}/scbsp.txt"
    shell:
        """
        Rscript scripts/scbsp.R --counts {input.counts} --coords {input.coords} --alpha {params.alpha} --output {output}
        """

rule smash:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv",
        smash_installed = "tools/SMASH-package-main/.installed"
    output:
        'results/{dataset}/smash_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/smash.yaml'
    benchmark:
        "benchmarks/{dataset}/smash.txt"
    shell:
        """
        python scripts/smash.py \
            --alpha {params.alpha} \
            --output {output} \
            --coords {input.coords} \
            --counts {input.counts} \
            --path {workflow.basedir}/tools/SMASH-package-main
        """

rule spagft:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spagft_svgs.csv'
    params:
        alpha = alpha_in
    conda:
        'envs/spagft.yaml'
    benchmark:
        "benchmarks/{dataset}/spagft.txt"
    shell:
        """
        python scripts/spagft.py \
            --alpha {params.alpha} \
            --output {output} \
            --coords {input.coords} \
            --counts {input.counts}
        """

rule spatialde:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spatialde_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/spatialde.yaml'
    benchmark:
        "benchmarks/{dataset}/spatialde.txt"
    shell:
        """
        python scripts/spatialde.py \
            --alpha {params.alpha} \
            --output {output} \
            --coords {input.coords} \
            --counts {input.counts}
        """

rule spatialde2:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spatialde2_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/spatialde2.yaml'
    benchmark:
        "benchmarks/{dataset}/spatialde2.txt"
    shell:
        """
        python scripts/spatialde2.py \
            --alpha {params.alpha} \
            --output {output} \
            --coords {input.coords} \
            --counts {input.counts}
        """

rule spark:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spark_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/r-env.yaml'
    benchmark:
        "benchmarks/{dataset}/spark.txt"
    shell:
        """
        Rscript scripts/spark-new.R --counts {input.counts} --coords {input.coords} --alpha {params.alpha} --output {output}
        """

rule sparkx:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/sparkx_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/r-env.yaml'
    benchmark:
        "benchmarks/{dataset}/sparkx.txt"
    shell:
        """
        Rscript scripts/sparkx.R --counts {input.counts} --coords {input.coords} --alpha {params.alpha} --output {output}
        """

rule spagene:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/spagene_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/r-env.yaml'
    benchmark:
        "benchmarks/{dataset}/spagene.txt"
    shell:
        """
        Rscript scripts/spagene.R --counts {input.counts} --coords {input.coords} --alpha {params.alpha} --output {output}
        """

rule meringue:
    input:
        counts = "data/{dataset}/counts.csv",
        coords = "data/{dataset}/coords.csv"
    output:
        'results/{dataset}/meringue_svgs.csv'
    params:
        alpha = alpha_in,
        dataset = "{dataset}"
    conda:
        'envs/r-env.yaml'
    benchmark:
        "benchmarks/{dataset}/meringue.txt"
    shell:
        """ 
        Rscript scripts/meringue.R --coords {input.coords} --counts {input.counts} --alpha {params.alpha} --output {output}
        """
