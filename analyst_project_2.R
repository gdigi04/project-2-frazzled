library(dplyr)

# temp_file_path = '/project/bf528/project_2/data/P4_vs_P7_cuffdiff_out/gene_exp.diff'

#read in cuffdiff output 
file_path = '/projectnb/bf528/users/frazzled/project_2/project-2-frazzled/cuffdiff_out/gene_exp.diff'	   
data = read.csv(file_path, sep = '')
ordered_data <- data[order(data$q_value),]
print(head(ordered_data))

# Produce a table of the top ten differentially expressed genes, 
# with their names, FPKM values, log fold change, p-value, and q-value in your report.
top_DEG = head(ordered_data, 10)
columns = c('gene_id', 'gene', 'value_1', 'value_2', 'log2.fold_change.', 'p_value', 'q_value')
report_table <- subset(top_DEG, select=columns)
write.csv(report_table, file = '/projectnb2/bf528/users/frazzled/project_2/top_DEG_table.csv')

#histogram of all log2 values
log2 = data$log2.fold_change.
hist(log2, breaks = 50)

#filter data to remove insignificant entries
sig_data = filter(data, significant == "yes") #threshold = 0.0075
hist(sig_data$log2.fold_change., breaks = 50)

#creating dataframes of upregulated and downregulated genes, respectively
upregulated = filter(sig_data, log2.fold_change. > 0)
downregulated = filter(sig_data, log2.fold_change. < 0)

#number of up & downregulated genes
up_num = nrow(upregulated)
down_num = nrow(downregulated)
print(up_num)
print(down_num)

write(upregulated$gene, file = '/projectnb2/bf528/users/frazzled/project_2/upregulated_genes.txt')
write(downregulated$gene, file = '/projectnb2/bf528/users/frazzled/project_2/downregulated_genes.txt')

downDAVID = read.csv(file = '/projectnb2/bf528/users/frazzled/project_2/downregulated_genes_DAVID.csv')
upDAVID = read.csv(file = '/projectnb2/bf528/users/frazzled/project_2/upregulated_genes_DAVID.csv')


