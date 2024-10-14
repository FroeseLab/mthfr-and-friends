# mthfr-and-friends
A repository for the ap ms analysis of MTHFR

Running order for scripts in folder "scripts":
1) "Script_SAINT_list_input_v1.R"
2) "Script_SAINT_output_analysis_v1.R"
3) "Script_Venn_Graph_v1.qmd"
4) "Script_SAINT_output_analysis_38-656_as_control.R"

SAINT analysis performed via online server https://reprint-apms.org/?q=analysis_front_apms
Input files found in folder "data/interim"
Main analysis using empty vector as control: "SAINT_list_input.csv"
Secondary analysis with MTHFR_38-656 as control: "SAINT_list_input_MTHFR38to656_control.csv"

Crapome retrieved from online server https://reprint-apms.org/?q=chooseworkflow 
Input  Protein Identifiers found in file "crapome_input.csv" in folder "data/interim"
