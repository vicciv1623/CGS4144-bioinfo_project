library(topGO)
library(org.Hs.eg.db)

res <- read.delim("results-daniel/step3_DE_results_all.tsv")

all_genes <- res$Gene 
sig_genes <- res$Gene[res$significant != "Not significant"]

gene_list <- factor(as.integer(all_genes %in% sig_genes))
names(gene_list) <- all_genes
print(table(gene_list))


#build the topGO object

GOdata <- new("topGOdata",
ontology = "BP",
allGenes = gene_list,
nodeSize = 10,
annot = annFUN.org,
mapping = "org.Hs.eg.db",
ID = "ensembl")

GOdata

result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
result

#collect results sorter=d= by p value
n_terms <- length(score(result))
tab <- GenTable(GOdata, pvalue = result, orderBy = "pvalue", topNodes = n_terms, numChar = 1000)

tab$pvalue <- as.numeric(sub("< ", "", tab$pvalue))

print(head(tab, 15))
cat("terms tested:", nrow(tab), "\n")
cat("terms with p < 0.01:", sum(tab$pvalue < 0.01), "\n")
write.csv(tab, "results-daniel/step5_topGo_BP.csv", row.names = FALSE)
