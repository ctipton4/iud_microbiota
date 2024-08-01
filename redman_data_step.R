library(ape)
library(phyloseq)
library(biomformat)
library(phytools)

## this code reads in the .biom file to start building a phyloseq object
# the biom file contains count information and taxonomic info at every level
db = read_biom("FullTaxa.otu_table.biom")
otu = otu_table(as.matrix(biom_data((db))), taxa_are_rows = T)
om = observation_metadata(db)$taxonomy
library(stringr)
om = str_split_fixed(om, " ; ", 7)
tt = tax_table(om)
rownames(tt) = 0:(dim(tt)[1]-1)
d = merge_phyloseq(otu, tt)
rm(db)
rm(tt)
rm(om)
rm(otu)

# ### these are custom functions you will probably need at some point 
# scripts = list.files("Custom_R_functions/")
# for (i in scripts){
#   source(file.path(paste("Custom_R_functions/", i, sep = "")))
# }

#this function is a custom script from
#http://www.phytools.org/read.newick/v0.5/read.newick.R
tre = read.newick("Redman_6469B.otus.tre")


#the tree and otu table hav different otu labels... they are linked in our mapping file
f = "Redman_6469B.otu_map.condensed.txt"
no_col <- max(count.fields(f, sep = "\t"))
map <- read.table(f,sep="\t",fill=TRUE,col.names=1:no_col)

names(map) = c("label","code")
map$code = gsub(":", "_", map$code)
map$code = gsub(";", "_", map$code)

l1 = data.frame(code = tre$tip.label)
lnew = merge(l1, map, by = "code", sort = F)

if(all(tre$tip.label == lnew$code))
   tre$tip.label = lnew$label

sample_names(d) = tolower(sample_names(d))
sample_names(d) = gsub('-', '.', sample_names(d))
sample_names(d) = gsub('_', '.', sample_names(d))

#midpoint rooting the tre
library(phytools)
tre = midpoint.root(tre)


#set up the metadata
met = data.frame(sample_names(d))
names(met)[1] = "sample"
row.names(met) = met$sample
met$sample = tolower(met$sample)

## custom metadata setup ## 
# create 'Patient' to link metadata and sequencing data
met$Patient = gsub(".357wf.806r", "", met$sample)

# read in and cleanup metadata if necessary
met1 = read.csv(file = "metadata.csv", header = T)

met1$Patient = gsub("IUD ", "", met1$ID)
dim(met)
dim(met1)
met = merge(met, met1, by = "Patient")
dim(met) # always check that dimensions make sense before and after mergers

## since Mirena is better represented, create a variable for comparing mirena to others
met$Mirena = met$Brand
met$Mirena = as.character(met$Mirena)
met$Mirena[which(met$Mirena != "Mirena")] = "Other"

## Calculate duration that IUD was used
met$Date.of.procedure = as.Date(as.character(met$Date.of.procedure), format = "%m/%d/%Y")
met$date.IUD.placed = as.Date(as.character(met$date.IUD.placed), format = "%m/%d/%Y")
met$Duration = as.numeric(met$Date.of.procedure - met$date.IUD.placed)
met$DurationWeeks = round(met$Duration/7, digits = 1) # change from days to months
met$DurationYears = round(met$Duration/365, digits = 1) # change from days to months

## drop variables that won't be used from here on out
# met$ID = met$Date.of.procedure = met$Ht = met$Wt = met$Gravida = met$Parity = met$date.IUD.placed = met$Date.samples.shipped = NULL
met$Date.of.procedure = met$Ht = met$Wt = met$date.IUD.placed = met$Date.samples.shipped = NULL

## use only hispanic
met = met[which(met$Ethnicity == "Hispanic"),]

## merge the OTU, tree, and metadata
row.names(met) = met$sample
d = merge_phyloseq(d, tre, sample_data(met))

#give the otu table the correct header
colnames(tax_table(d)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#just to help distinguish taxa ids from numbers
taxa_names(d) = paste("OTUid.", taxa_names(d), sep = '')

#remove informationless OTU entries
d <- prune_taxa(taxa_sums(d) > 10, d) # drop rare OTUs receiving less than 10 reads
d  = subset_taxa(d, Kingdom != "No Hit")
d  = subset_taxa(d, Genus != "Ralstonia") # common reagent contaminant
# d  = subset_taxa(d, Genus != "Halospirulina")
# d  = subset_taxa(d, Genus != "Methylobacterium")
# d  = subset_taxa(d, Genus != "Stenotrophomonas")
# d  = subset_taxa(d, Genus != "Sphingomonas")
# d  = subset_taxa(d, Genus != "Elizabethkingia")
# d  = subset_taxa(d, Family != "Sphingomonadaceae")
# d  = subset_taxa(d, Order != "Oscillatoriales")
# d  = subset_taxa(d, Order != "Burkholderiales")
# d  = subset_taxa(d, Genus != "Bacillus")

#using this to assess read depth, and removing those that are less that 500 reads (454 standard)
x = data.frame(t(otu_table(d)))
temp = data.frame(rowSums(x)) # all look good
temp$fullname = rownames(temp)
# list = temp[temp$rowSums.x. > 500,]
# list = row.names(list)
# d = prune_samples(list, d)


#sample names can have no weird characters and must start with a letter
######calc gen_p and unif data

otu = data.frame(t(otu_table(d)))

# #create trimmed family file and pd
# d_s <- tax_glom(d, taxrank = "Family")
# pd = sample_data(d_s)
# #dp = transform_sample_counts(dp, function(x) x/sum(x))
# gc = data.frame(t(otu_table(d_s)))
# 
# n = paste(tax_table(d_s)[,1]," ; ", 
#           tax_table(d_s)[,2]," ; ",
#           tax_table(d_s)[,3]," ; ",
#           tax_table(d_s)[,4]," ; ",
#           tax_table(d_s)[,5]," ; ",
#           tax_table(d_s)[,6], sep = '')
# n = clean_blast2_names(n, tax_level="family")
# 
# names(gc) = n
# fam_p = gc/rowSums(gc)
# 
# #create trimmed phylum file and pd
# d_p <- tax_glom(d, taxrank = "Phylum")
# pd = sample_data(d_p)
# #dp = transform_sample_counts(dp, function(x) x/sum(x))
# gc = data.frame(t(otu_table(d_p)))
# 
# n = paste(tax_table(d_p)[,1]," ; ", 
#           tax_table(d_p)[,2]," ; ",
#           tax_table(d_p)[,3]," ; ",
#           tax_table(d_p)[,4]," ; ",
#           tax_table(d_p)[,5]," ; ",
#           tax_table(d_p)[,6], sep = '')
# n = clean_blast2_names(n, tax_level="phylum")
# 
# names(gc) = n
# phy_p = gc/rowSums(gc)

#create trimmed genus file and pd
d_s <- tax_glom(d, taxrank = "Genus")
pd = sample_data(d_s)
#dp = transform_sample_counts(dp, function(x) x/sum(x))
gc = data.frame(t(otu_table(d_s)))

n = paste(tax_table(d_s)[,1]," ; ", 
          tax_table(d_s)[,2]," ; ",
          tax_table(d_s)[,3]," ; ",
          tax_table(d_s)[,4]," ; ",
          tax_table(d_s)[,5]," ; ",
          tax_table(d_s)[,6], sep = '')
n = clean_blast2_names(n, tax_level="genus")

names(gc) = n
gen_p = gc/rowSums(gc)



#create trimmed species file and pd
d_s <- tax_glom(d, taxrank = "Species")
pd = sample_data(d_s)
#dp = transform_sample_counts(dp, function(x) x/sum(x))
gc = data.frame(t(otu_table(d_s)))

n = paste(tax_table(d_s)[,1]," ; ", 
          tax_table(d_s)[,2]," ; ",
          tax_table(d_s)[,3]," ; ",
          tax_table(d_s)[,4]," ; ",
          tax_table(d_s)[,5]," ; ",
          tax_table(d_s)[,6]," ; ",
          tax_table(d_s)[,7], sep = '')
n = clean_blast2_names(n, tax_level="species")

names(gc) = n
spp_p = gc/rowSums(gc)

# library(doParallel)
# cl <- makeCluster(3, rscript_args = c("--no-init-file", "--no-site-file","--no-environ"))
# #for linux   .... cl <- makePSOCKcluster(2)
# registerDoParallel(cl)
unif_w = UniFrac(d, weighted=T, parallel=F, normalized = T)
unif_uw = UniFrac(d, weighted=F, parallel=F)
# stopCluster(cl)

unif_w = as.matrix(unif_w)
unif_uw = as.matrix(unif_uw)

