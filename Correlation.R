#correlation
cor_matrix <- cor(df[, c(2:11, 21, 29, 37, 45, 53, 61, 69, 77, 85, 93)], use = "pairwise.complete.obs")

library(ggplot2)
library(reshape2)

melted_cormat <- melt(cor_matrix)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) + 
  coord_fixed()
##FINISHED FOR HYBRID##

##NOW LETS DO IT FOR EVERY OTHER METHOD##

# Sample attributes separated for Illumina and ONT
sample_attributes_ilu <- df[, 2:6]
sample_attributes_ont <- df[, 7:11]

# Function to extract method columns for ILU methods
get_method_data_ilu <- function(start_col) {
  method_columns <- seq(start_col, by=8, length.out=10)
  return(cbind(sample_attributes_ilu, df[, method_columns]))
}

# Function to extract method columns for ONT methods
get_method_data_ont <- function(start_col) {
  method_columns <- seq(start_col, by=8, length.out=10)
  return(cbind(sample_attributes_ont, df[, method_columns]))
}

ilu_ori_sph <- get_method_data_ilu(14)
ilu_alt_RKI <- get_method_data_ilu(15)
ilu_alt2_man <- get_method_data_ilu(16)
ont_ori_mbi <- get_method_data_ont(17)
ont_alt_ART <- get_method_data_ont(18)
ont_alt2_man <- get_method_data_ont(19)
ilu_dnv_spa <- get_method_data_ilu(20)
hybrid <- cbind(sample_attributes_ilu, sample_attributes_ont, df[, seq(21, by=8, length.out=10)])

cor_ilu_ori_sph <- cor(ilu_ori_sph, use = "pairwise.complete.obs")
cor_ilu_alt_RKI <- cor(ilu_alt_RKI, use = "pairwise.complete.obs")
cor_ilu_alt2_man <- cor(ilu_alt2_man, use = "pairwise.complete.obs")
cor_ont_ori_mbi <- cor(ont_ori_mbi, use = "pairwise.complete.obs")
cor_ont_alt_ART <- cor(ont_alt_ART, use = "pairwise.complete.obs")
cor_ont_alt2_man <- cor(ont_alt2_man, use = "pairwise.complete.obs")
cor_ilu_dnv_spa <- cor(ilu_dnv_spa, use = "pairwise.complete.obs")
cor_hybrid <- cor(hybrid, use = "pairwise.complete.obs")

visualize_cor_matrix(cor_ilu_ori_sph, "ilu_ori_sph")
visualize_cor_matrix(cor_ilu_alt_RKI, "ilu_alt_RKI")
visualize_cor_matrix(cor_ilu_alt2_man, "ilu_alt2_man")
visualize_cor_matrix(cor_ont_ori_mbi, "ont_ori_mbi")
visualize_cor_matrix(cor_ont_alt_ART, "ont_alt_ART")
visualize_cor_matrix(cor_ont_alt2_man, "ont_alt2_man")
visualize_cor_matrix(cor_ilu_dnv_spa, "ilu_dnv_spa")
visualize_cor_matrix(cor_hybrid, "hybrid")

install.packages("heatmaply")
library(heatmaply)

heatmaply_cor(cor_ilu_ori_sph, main = "ILU ORI SPH")
heatmaply_cor(cor_hybrid, main = "HYBRID", colors = )

cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(mtcars)
