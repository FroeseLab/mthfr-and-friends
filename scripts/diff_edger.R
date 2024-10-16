library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(stringr)


# Load the data
dat_raw <- read.csv("data/raw/Prepare_SAINT_list_input_file/SAINT_list_input.csv")

data_raw_scaffold<- read.csv("data/raw/Prepare_SAINT_list_input_file/Proteins_Report_raw.csv", skip = 1, header = TRUE, check.names = FALSE)

# Get metadta from the raw data
dat_raw %<>%
    mutate(replicate = str_match( X.APName, '[0-9+]F$')[,1]) %>%
    mutate(gene=translate_up(X.PreyName)) %>%
    mutate(condition=X.BaitName_Condition)
# Reformat the data into a count matrix
counts = dat_raw %>%
    mutate(bait_name = gsub('_','',X.BaitName_Condition)) %>%
    mutate(sample_name = paste(bait_name, X.APName, sep="_")) %>%
    tidyr::pivot_wider(id_cols = X.PreyName,
                       names_from = sample_name,
                       values_from = X.SpectralCounts,
                       values_fill=0) %>%
    tibble::column_to_rownames(var = "X.PreyName")

# Extract experimental groups
groups = factor(sapply(strsplit(colnames(counts), "_"), function(x) x[1]),
                levels = c("CONTROL" ,  "MTHFR" ,"MTHFRT34A" ,  "MTHFRtrunc"  ))



# generate a DGEList object
y = DGEList(counts, groups=groups)

# Estimate dispersion: note that 
# the dispersion estimates will be overwritten later
# as for AP-MS I think we should have a constant library size
# for all the samples
y = estimateDisp(y)

# plot MDS - Multidimensional scaling plot
plotMDS(y)

# Overwrite the library size with the mean library of all samples
y$samples$lib.size <- mean(y$samples$lib.size)

# Create a design matrix from the groups
design <- model.matrix(~groups)
# Fit the model
fit <- glmQLFit(y, design)

# Compare MTHFR vs CONTROL
qlf.2vs1 <- glmQLFTest(fit, coef=2)
head(qlf.2vs1$table)

ggplot(qlf.2vs1$table) +
    geom_point(
            aes(x=logFC,
            y=-log10(PValue),
            color=qvalue(p=PValue)$qvalues < 0.1 & abs(logFC) > 5),
    ) +
    ggtitle("MTHFR vs CONTROL")

# Compare MTHFR vs truncate
mthfr_vs_trunc <- glmQLFTest(fit, contrast=c(0,-1,0,1))

ggplot(mthfr_vs_trunc$table) +
    geom_point(
            aes(x=logFC,
            y=-log10(PValue),
            color=qvalue(p=PValue)$qvalues < 0.1 & abs(logFC) > 5),
    ) +
    ggtitle("MTHFR vs truncated does not show significant changes")


# Test all vs control
vs_ctrl <- glmQLFTest(fit, coef=2:4)


vs_ctrl = add_genes(vs_ctrl)
hits = topTags(vs_ctrl, sort.by="PValue", n=400)$table

# Lots of hits
ggplot(vs_ctrl$table) +
    geom_point(
            aes(x=logFC.groupsMTHFR,
            y=-log10(PValue),
            color=qvalue(p=PValue)$qvalues < 0.1 & abs(logFC.groupsMTHFR) > 5),
    ) +
    ggtitle("CONTROL vs all")


# Check raw data

dat_raw %>%
    filter(gene %in% c( 'MTHFR')) %>%
    
    ggplot(aes(y=X.SpectralCounts,
        x=X.BaitName_Condition,
        color=replicate)) +
    facet_wrap(~gene) +
    geom_point() +
    scale_y_log10()


dat_raw %>%
    filter(gene %in% c( 'PRKDC')) %>%
    
    ggplot(aes(y=X.SpectralCounts,
        x=X.BaitName_Condition,
        color=replicate)) +
    facet_wrap(~gene) +
    geom_point() +
    scale_y_log10()


# Check number of genes in control
dat_raw %>%
    filter(condition == 'CONTROL') %>%
    pull(gene) %>%
    unique() %>%
    length()

    
dat_raw %>%
    merge(head(hits, 132)) %>%
    mutate(p_val = forcats::fct_rev(format(log10(PValue), digits=3))) %>%
    ggplot(aes(y=X.SpectralCounts,
        x=X.BaitName_Condition,
        color=replicate)) +
    facet_wrap(~ p_val + gene) +
    geom_point() +
    scale_y_log10() +
    ggtitle("Top 132 hits") +
    # rotate axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))



dat_raw %>%
    filter(condition == 'CONTROL') %>%
     ggplot(aes(y=X.SpectralCounts,
        x=gene,
        #color=replicate,
        color =  gene == 'MTHFR')) +
   
    geom_point() +
    ggtitle('MTHFR epressed in control') +
    scale_y_log10()



# Directly contrast MTHFR vs truncated
mthfr_vs_trunc = add_genes(mthfr_vs_trunc)
topTags(mthfr_vs_trunc, sort.by="PValue", n=20)

# Directly contrast MTHFR vs mutated
mthfr_vs_mut <- glmQLFTest(fit, contrast=c(0,-1,1,0))
mthfr_vs_mut = add_genes(mthfr_vs_mut)
topTags(mthfr_vs_mut, sort.by="PValue", n=20)

# -> No significance for both
