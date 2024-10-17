rule all:
    input:
        "data/output/saint_output_merged_main.csv",


rule prepare_saint:
    conda:
        "environment.yml"
    input:
        data_raw_scaffold="data/raw/Prepare_SAINT_list_input_file/Proteins Report of merged human files with all 4 conditions (1. Feb, 5. Feb, 18. Feb).txt",
        data_rename="data/metadata/SAINT_rename_file.csv",
        renv="restored.Rprofile",
    output:
        data_saint="data/interim/SAINT_list_input.csv",
        data_saint_MTHFR38to65="interim/SAINT_list_input_MTHFR38to656_control.csv",
        crapome="data/interim/crapome_input.csv",
    script:
        "scripts/Script_SAINT_list_input_v1.R"


rule run_saint:
    message:
        """
        Data from SAINT analysis using "data/interim/SAINT_list_input.csv" as input, where empty vector is control
        SAINT analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server 
        Import raw .tsv file of the "raw SAINT results".
        When downloaded a .output file is given.
        Just add .tsv file format

        Bait-prey output was exported from the unfiltered matrix
        from "Bait-prey heatmap" tab. Downloaded as a .tab file.
        Just add .tsv file format
        """
    input:
        data_saint="data/interim/SAINT_list_input.csv",
        saint_output="data/resources/prepared/21262.output.tsv",
        bait_prey_output="data/resources/prepared/21262_baitPrey_data.tab.tsv",
    output:
        saint_output="data/interim/saint_output_main.csv",
        bait_prey_output="data/interim/saint_baitPrey_main.tab.tsv",
    shell:
        """
        cp {input.saint_output} {output.saint_output}
        cp {input.bait_prey_output} {output.bait_prey_output}
        """


rule run_saint_38656:
    message:
        """
        Data from SAINT analysis using "data/interim/SAINT_list_input_MTHFR38to656_control.csv" as input, where empty vector is control
        SAINT analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server 
        Import raw .tsv file of the "raw SAINT results".
        When downloaded a .output file is given.
        Just add .tsv file format

        Bait-prey output was exported from the unfiltered matrix
        from "Bait-prey heatmap" tab. Downloaded as a .tab file.
        Just add .tsv file format
        """
    input:
        data_saint="data/interim/SAINT_list_input_MTHFR38to656_control.csv",
        saint_output="data/resources/prepared/21302.output.tsv",
        bait_prey_output="data/resources/prepared/21302_baitPrey_data.tab.tsv",
    output:
        saint_output="data/interim/saint_output_MTHFR38to656.csv",
        bait_prey_output="data/interim/saint_baitPrey_MTHFR38to656.tab.tsv",
    shell:
        """
        cp {input.saint_output} {output.saint_output}
        cp {input.bait_prey_output} {output.bait_prey_output}
        """


rule merge_saint_outputs:
    input:
        saint_output="data/interim/saint_output_{condition}.csv",
        saint_prey_output="data/interim/saint_output_{condition}.csv",
    output:
        saint_output="data/output/saint_output_merged_{condition}.csv",
    script:
        "scripts/merge_saint_outputs.R"


# Use renv
wildcard_constraints:
    # Allows renv_dir also to be the empty (=root dir)
    renv_dir=".*",


rule restore_renv:
    conda:
        "environment.yml"
    input:
        rprofile="{renv_dir}.Rprofile",
        lock="{renv_dir}renv.lock",
    output:
        restored="{renv_dir}restored.Rprofile",
    log:
        "logs/renv/{renv_dir}_restore.log",
    shell:
        """
        R -e 'source("{input.rprofile}", chdir=TRUE);renv::restore()'
        echo 'source("{input.rprofile}", chdir=TRUE)' > {output.restored}
        """
