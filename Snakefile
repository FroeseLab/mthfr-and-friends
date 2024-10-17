rule prepare_saint:
    conda:
        "environment.yml"
    input:
        data_raw_scaffold =  "data/raw/Prepare_SAINT_list_input_file/Proteins Report of merged human files with all 4 conditions (1. Feb, 5. Feb, 18. Feb).txt",
        data_rename = "data/metadata/SAINT_rename_file.csv"
    output:
        data_saint = "data/interim/SAINT_list_input.csv",
        data_saint_MTHFR38to65 = "interim/SAINT_list_input_MTHFR38to656_control.csv",
        crapome = "data/interim/crapome_input.csv"
    script:
        "scripts/Script_SAINT_list_input_v1.R"
