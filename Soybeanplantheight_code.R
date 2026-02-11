###GWAS for soybean plant height 

install.packages("rMVP")
library(rMVP)
#data preparation and input
MVP.Data(fileVCF="SNPs_onlyfiltered_soybean.vcf",
         filePhe="StorageID(GWAS).txt",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.vcf")
genotype <- attach.big.matrix("mvp.vcf.geno.desc")
phenotype <- read.table("mvp.vcf.phe",head=TRUE)
map <- read.table("mvp.vcf.geno.map" , head = TRUE)
##adding environment
Envt_data <- read.csv("Envt_data.csv")
str(Envt_data)
Covariates <- model.matrix.lm(~as.factor(method_name), data=Envt_data, na.action = "na.pass")
####GWAS
imMVP <- MVP(
  phe=phenotype,
  geno=genotype,
  map=map,
  CV.GLM=Covariates,     
  CV.MLM=Covariates,
  CV.FarmCPU=Covariates,
  nPC.GLM=5,      
  nPC.MLM=3,  
  nPC.FarmCPU=3,
  maxLine=100,          
  ncpus=4,
  vc.method="BRENT",      
  method.bin="static",  
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

################################################################################################
##Machine learning for soybean plant height

#phase impute using beagle
java -jar /software/projects/pawsey0149/mnatukunda/setonix/beagle/beagle.28Jun21.220.jar gt=prunedSNPs.vcf.gz out=phased_SNPs impute=true 
#encoding genotype data using vcf tools
vcftools --gzvcf PhasedSNPs.vcf.gz --012 --out PhasedSNPs_encoded
mle <- read_csv("EncodedML.csv")
pheno <- read_csv("StorageID_GWAS.csv")
rownames(mle) <- pheno$Ind
#add plant height and environment
mle$pheno <- pheno$Pheno 
envt <- read_csv("Method_data.csv") 
envt1 <- envt %>% 
  select(-Data.storage.ID)
str(envt1)
envt1$method_name <- as.factor(envt1$method_name)
str(envt1)
dim(envt1)
#convert environment data to numeric format 
en_matrix <- model.matrix( ~ . -1, data = envt1)
dim(en_matrix)
#merge environment and encoded data 
merged_data <- bind_cols(mle, en_matrix)
dim(merged_data)
##split data into 80:20 ratio
set.seed(40) 
datasize <- floor(0.8*nrow(mle)) 
train_indices <- sample(seq_len(nrow(merged_data)), size = datasize) 
train_data <- merged_data[train_indices, ]
dim(train_data) 
train1 <- data.frame(train_indices)
write_csv(train1, file = "train_indices", col_names = TRUE)
write_csv(train_data, file = "train_data.csv", col_names = TRUE)
validation_data <- merged_data[-train_in, ] 
write_csv(validation_data, file = "validation_data.csv", col_names = TRUE)
dim(validation_data) 
head(validation_data)
##separate features from response variables
train_data <- read_csv("train_data.csv")
featurest <- train_data %>%
  select(-pheno)
str(featurest)
labelt <- train_data$pheno
is.data.frame(featurest)
is.data.frame(labelt)
str(labelt)
#train model with xgboost
library(xgboost)
featurest_matrix <- as.matrix(featurest)
dim(featurest_matrix)

set.seed(30) 
##parameter tuning using cross validation and loop
generate_params <- function() {
  list(
    max_depth = sample(3:10, 1),
    eta = runif(1, 0.01,0.3),
    gamma = runif(1,0,5),
    subsample = runif ( 1, 0.6, 1),
    colsample_bytree = runif(1,0.5, 1),
    min_child_weight = sample(1:10, 1),
    objective = "reg:squarederror",
    eval_metric = "rmse" ) }
results <- list()
n_iter <- 10

for (i in 1:n_iter) {
  params <- generate_params()
  
  cv_results <- xgb.cv(params = params, data = dtrain, nrounds = 100, nfold= 3, early_stopping_rounds = 10,
                       verbose = 0)
  
  results[[i]] <- list(
    params = params,
    best_iteration = cv_results$best_iteration,
    best_error = min(cv_results$evaluation_log$test_error_mean)) }

object.size(dall)
length(unlist(folds[-k]))
##combine results
library(data.table)
results_dt <- rbindlist(lapply(results, function(x) {
  c(x$params, best_iteration = x$best_iteration, best_error = x$best_error)
}), fill = TRUE)
##best parameters
best_pm <- results_dt[order(best_error)][1]
print(best_pm)
#train with optimum number of rounds and best parameters
set.seed(100)
bst <- xgboost(data = as.matrix(featurest), 
               label = labelt, 
               booster = "gbtree", 
               objective = "reg:squarederror", 
               nrounds = 28, 
               verbose = 0, 
               eta = 0.189, 
               max_depth = 3, 
               subsample = 0.763, 
               colsample_bytree = 0.926, 
               reg_alpha = 0, 
               reg_lambda = 1, 
               gamma = 1.967, 
               min_child_weight = 2)

#check model performance
validation_data <- read_csv("validation_data.csv")
dim(validation_data)
valph <- validation_data$pheno
vall2 <- validation_data %>%
  select(-pheno)
###performance metrics on training data set  
library(caret)
preds1 <- predict(bst, as.matrix(featurest), type = "response")  
postResample(pred = preds1, obs = labelt) 
# performance metrics on the test/validation data 
preds <- predict(bst, as.matrix(vall2), type = "response")  
summary(valph) 
postResample(pred = preds, obs = valph) 

##feature importance to obtain top influencing SNPs
importance_matrix <- xgb.importance(model = bst) 
print(importance_matrix) 
head(importance_matrix, n = 25) 
ftimp <- head(importance_matrix, n = 20) 
write_csv(ftimp, file = "top20SNPs_ML.csv", col_names = TRUE)

xgb.plot.importance(importance_matrix = importance_matrix, top_n = 10, xlab = 'Gain' )  

xgb.ggplot.importance(importance_matrix = importance_matrix, top_n = 10) +
  ggtitle("Top 10 Feature Importance") + 
  theme_minimal()

#############################################################################################
##Local haplotyping for soybean plant height
install.packages("devtools")
devtools::install_github("jacobimarsh/crosshap")
library(crosshap)
library(tidyverse)

#local haplotyping for chromosome 19 QTL (92KB region with 313 SNPs)
vcf <- read_vcf('92kbregion.vcf')
LD <- read_LD('92kbregionLDmatrix.ld', vcf = vcf)
pheno <- read_pheno ('StorageID(GWAS).txt')
metadata <- read_metadata('Origin_data.txt')
str(metadata)
metadata$Metadata <- dplyr::recode(metadata$Metadata, 
                                   Modern_cultivar = "Modern cultivar",
                                   Old_cultivar = "Old cultivar",
                                   Wild_soybean = "Wild soybean")
head(vcf)
head(LD, c(4,4))
head(pheno)
head(metadata)
MGmin <- 30 
epsilon <- c(0.2,0.4,0.6,0.8,1) 
HapObject <- run_haplotyping(vcf = vcf,
                             LD = LD,
                             pheno = pheno,
                             metadata = metadata,
                             epsilon = epsilon,
                             MGmin = MGmin) 
hap_clustree <- clustree_viz(HapObject,
                             type = 'hap')
hap_clustree
#visualizing
Hap_viz <- crosshap_viz(HapObject = HapObject, epsilon = 0.6)
Hap_viz
#swap hap results obtained to get the chromosomal position of each SNP
Hap_viz <- crosshap_viz(HapObject = HapObject, epsilon = 0.6, plot_left = 'pos')
Hap_viz
hapfile <- HapObject$Haplotypes_MGmin30_E0.6$Hapfile
hapfile

#reorganizing individual plots
library(patchwork)
top <- build_top_metaplot(HapObject = HapObject, epsilon = 0.6, hide_labels = F)
top
bot <- build_bot_halfeyeplot(HapObject = HapObject, epsilon = 0.6, hide_labels = T,)
bot
left <- build_left_posplot(HapObject = HapObject, epsilon = 0.6, hide_labels = F)
left
middle <- build_mid_dotplot(HapObject = HapObject, epsilon = 0.6, hide_labels = F)
middle
rightphenoplot<- build_right_phenoplot(HapObject = HapObject, epsilon = 0.6, hide_labels = F)
rightphenoplot
#Extract and build summary tables as images 
library(ggplotify)
summaryMGtables <- build_summary_tables(HapObject = HapObject, epsilon= 0.6 )
MGsummarytable <- summaryMGtables[[1]]
Hapssummarytable <- summaryMGtables[[2]]
MGsplot <- as.ggplot(MGsummarytable) +
  theme_void()+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),plot.background = element_blank(),
        panel.background = element_blank())+
  coord_cartesian(expand = TRUE)

hapsplot <- as.ggplot(Hapssummarytable) +
  theme_void() +  # Remove titles or axis labels
  theme(plot.margin = margin(t = 0, r = 60, b = 0, l = 0))  # Add right margin
hapsplot
#Wrapping plots together with patchwork
final_plot10_3 <- (
  plot_spacer() +                    # Spacer (optional)
    top +                           # Top 20 plot
    wrap_elements(MGsplot) +         # MGtable
    left +                    # Left pheno plot
    middle +                     # Dot plot
    rightphenoplot +                   # Right pheno plot
    guide_area() +                     # Guide area
    bot +                           # Bottom 20 plot
    wrap_elements(hapsplot) +        # Haptable
    plot_layout(ncol = 3, guides = 'collect') +  # Keep the grid size
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &  # Custom tags
    theme(plot.tag = element_text(face = 'bold', size = 16), 
          legend.text = element_text(size= 8.5))  # Define grid layout
)


final_plot10_3
#saving plot 
ggsave(filename = "5crosshap_figure_newchr19.png", plot = final_plot10_3, width = 16, height = 10, units = "in", dpi = 600)

#extracting individual haplotype information 
chr19_93kbhaplotypes <- HapObject$Haplotypes_MGmin30_E0.6$Indfile
write.csv(chr19_93kbhaplotypes, file = "chr19_93kbMG_haplotypes", quote = F)
markergroups <- HapObject$Haplotypes_MGmin30_E0.6$Hapfile
mG_SNPs<- HapObject$Haplotypes_MGmin30_E0.6$Varfile 
write.csv(mG_SNPs, file = "chr19_93kbMG_SNPs", quote = F)

#calculating p.value threshold for A, B,D vs C
chr19_93kbhaplotypes$hap <- as.factor(chr19_93kbhaplotypes$hap)
chr19_93kbhaplotypes$Metadata <- as.factor(chr19_93kbhaplotypes$Metadata)
str(chr19_93kbhaplotypes)
boxplot(Pheno~chr19_93kbhaplotypes$hap, data = chr19_93kbhaplotypes)
Ntassigned <- subset(chr19_93kbhaplotypes,hap == "0")
Ntassigned
haplotypesF<- chr19_93kbhaplotypes %>% 
  filter(hap != "0")
head(haplotypesF)
str(haplotypesF)
boxplot(Pheno~hap, data = haplotypesF)
data3sampletest <- haplotypesF %>% 
  mutate(comparison_group = case_when(
    hap %in% c("A","B","D")~ "A",
    hap %in% c("C")~"C",
  ))
str(data3sampletest)
data3sampletest$comparison_group <- as.factor(data3sampletest$comparison_group)
t.test(Pheno~comparison_group, var.equal = FALSE, data = data3sampletest)

######################################################################################
#local haplotyping at chromosome 14 QTL region (18.3kb)
library(crosshap)
library(ggdist)
vcf <- read_vcf('chr14_18.3kb.vcf')
LD <- read_LD('chr14_18.3kbLDmatrix.ld', vcf = vcf)
pheno <- read_pheno ('StorageID(GWAS).txt')
metadata <- read_metadata ('Origin_data.txt')
metadata$Metadata <- dplyr::recode(metadata$Metadata, 
                                   Modern_cultivar = "Modern cultivar",
                                   Old_cultivar = "Old cultivar",
                                   Wild_soybean = "Wild soybean")
str(metadata)
head(vcf)
head(LD, c(4,4))
head(pheno)
head(metadata)
MGmin <- 20
epsilon <- c(0.2,0.4,0.6,0.8,1)
HapObject1 <- run_haplotyping(vcf = vcf,
                              LD = LD,
                              pheno = pheno,
                              metadata = metadata,
                              epsilon = epsilon,
                              MGmin = MGmin)
hap_clustree <- clustree_viz(HapObject1,
                             type = 'hap')
hap_clustree
#visualizing
Hap_viz1 <- crosshap_viz(HapObject = HapObject1, epsilon = 0.6)
Hap_viz1
###swap hap results obtained to get the chromosomal position of each SNP
Hap_viz1 <- crosshap_viz(HapObject = HapObject1, epsilon = 0.6, plot_left = 'pos')
Hap_viz1
ggsave("chr1418.3kb1crosshap.png", plot = Hap_viz1, width = 5000,height = 3700, 
       units= "px")
##############extract haplotype information
HapObject$Haplotypes_MGmin20_E0.6$Indfile 
Hapfile <- HapObject1$Haplotypes_MGmin20_E0.6$Hapfile
Hapfile
print(HapObject1$Haplotypes_MGmin20_E0.6$Varfile, n= 50) 
chr14_18.3kbmG_SNPs<- HapObject$Haplotypes_MGmin20_E0.6$Varfile 
write.csv(chr14_18.3kbmG_SNPs, file = "chr14_18.3kbmG_SNPs", quote = F)

top <- build_top_metaplot(HapObject = HapObject1, epsilon = 0.6, hide_labels = F)
top
bot <- build_bot_halfeyeplot(HapObject = HapObject1, epsilon = 0.6, hide_labels = T,)
bot
left <- build_left_posplot(HapObject = HapObject1, epsilon = 0.6, hide_labels = F)
left
middle <- build_mid_dotplot(HapObject = HapObject1, epsilon = 0.6, hide_labels = F)
middle
rightphenoplot<- build_right_phenoplot(HapObject = HapObject1, epsilon = 0.6, hide_labels = F)
rightphenoplot
##Extract and build summary tables as images
summaryMGtables <- build_summary_tables(HapObject = HapObject1, epsilon= 0.6 )
MGsummarytable <- summaryMGtables[[1]]
Hapssummarytable <- summaryMGtables[[2]]
MGsplot <- as.ggplot(MGsummarytable) +
  theme_void()+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),plot.background = element_blank(),
        panel.background = element_blank())+
  coord_cartesian(expand = TRUE)

hapsplot <- as.ggplot(Hapssummarytable) +
  theme_void() +  # Remove titles or axis labels
  theme(plot.margin = margin(t = 0, r = 60, b = 0, l = 0))  # Add right margin
hapsplot
###Wrapping plots together with patchwork
final_plot10_3 <- (
  plot_spacer() +                    # Spacer (optional)
    top +                           # Top plot
    wrap_elements(MGsplot) +         # MGtable
    left +                    # Left pheno plot
    middle +                     # Dot plot
    rightphenoplot +                   # Right pheno plot
    guide_area() +                     # Guide area
    bot +                           # Bottom plot
    wrap_elements(hapsplot) +        # Haptable
    plot_layout(ncol = 3, guides = 'collect') +  # Keep the grid size
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &  # Custom tags
    theme(plot.tag = element_text(face = 'bold', size = 16), 
          legend.text = element_text(size= 8.5))  # Define grid layout
)
final_plot10_3
##saving plot 
ggsave(filename = "crosshap_figure_newchr14.png", plot = final_plot10_3, width = 16, height = 10, units = "in", dpi = 600)

#calculating two sample t.test for haplotypes A and C vs B
help("t.test")
chr14hapcombination <- HapObject$Haplotypes_MGmin20_E0.6$Indfile
str(chr14hapcombination)
chr14hapcombination$hap <- as.factor(chr14hapcombination$hap)
chr14hapcombination$Metadata <- as.factor(chr14hapcombination$Metadata)
str(chr14hapcombination)
boxplot(Pheno~chr14hapcombination$hap, data = chr14hapcombination)
Ntassigned <- subset(chr14hapcombination,hap == "0")
Ntassigned
haplotypesF <- chr14hapcombination %>% ####exclude individuals not assigned haplotypes
  filter(hap != "0")
head(haplotypesF)
str(haplotypesF)
boxplot(Pheno~hap, data = haplotypesF)
data3sampletest <- haplotypesF %>% 
  mutate(comparison_group = case_when(
    hap %in% c("A","C")~ "A",
    hap %in% c("B")~"B",
  ))
str(data3sampletest)
data3sampletest$comparison_group <- as.factor(data3sampletest$comparison_group)
t.test(Pheno~comparison_group, var.equal = FALSE, data = data3sampletest)