##################################
# Functional Enrichment Analysis #
# using topGO                    #
##################################

library("topGO")

####################
# Data preparation #
####################
# Import gene-to-GO mappings: cotton GO annotation file downloaded from phytozome.
geneID2GO <- readMappings(file = "cotton.v2.1.stable.txt")

# Import all Gorai IDs as "universe"; note that for most enrichment analysis, "universe" should be background genes (e.g., genes included for DE or network analysis) instead of All genome.
id<-read.table("Dgenome_id.txt",header=TRUE,colClasses="character")
goraiUniverse<-id[,1]

# Here, we use goraiUniverse as beackground genes
universe<- goraiUniverse
head(universe)

# Gene list of interest could be a collection of DE genes, or co-expression module members
# but for now, I randomly sample 100 genes  
genes = universe[sample(1:length(universe),100)]
head(genes)

# locate genes in the university
geneList <- factor(as.integer(universe %in% genes))
names(geneList) <- universe
str(geneList)

#####################
# Make topGO object #
#####################
## Make Biological Process GOData object ##
# ontology: character string specifying the ontology of interest
#           ('BP':Biological process,'MF':molecular function or 'CC':cellular components)
# allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The
#           genes listed in this object define the gene universe.
# annotationFun: function that maps gene identifiers to GO terms. annFUN.gene2GO extracts the gene-to-GO mappings from customized object: geneID2GO
# nodeSize: nodes with at least  how many  genes, default is 1

GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
## display the GOdata object
GOdata.MF


############################
# Run enrichment analysis  #
############################
# fisher test
result<- runTest(GOdata.MF, algorithm = "classic", statistic = "fisher")
# Kolmogorov-Smirnov test.
resultKS <- runTest(GOdata.MF, algorithm = "classic", statistic = "ks")
# Kolmogorov-Smirnov test.
resultKS.elim <- runTest(GOdata.MF, algorithm = "elim", statistic = "ks")

##############################
# other algorithms supported #
##############################
# other algorithms that consider the GO structure-'elim' and 'weight':
r.elim = runTest(GOdata.MF, algorithm = 'elim', statistic = 'fisher')
r.weight = runTest(GOdata.MF, algorithm = 'weight', statistic = 'fisher')

# When the input list is compared with the previously computed background, or is a subset of reference list, choose hypergeometric or fisher, for latter only when your query number is quite small. When the input list has few or no intersections with the reference list, the Chi-square tests are more appropriate

whichAlgorithms()
# [1] "classic"     "elim"        "weight"      "weight01"    "lea"         "parentchild"

whichTests()
# [1] "fisher"     "ks"         "t"          "globaltest" "sum"        "ks.ties"


# get the enrichment p-values from all three algorithms and compare
pValue.classic <- score(r1.BP)
pValue.elim <- score(r1.BP.elim)[names(pValue.classic)]
pValue.weight <- score(r1.BP.weight)[names(pValue.classic)]
# for plots of p-values: set plotting parameters:
gstat <- termStat(GOdata.BP, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# function for building color map for plotting
colMap <- function(x) {
    .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
    return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol,log='xy')
plot(pValue.classic, pValue.weight, xlab = "p-value classic", ylab = "p-value weight",pch = 19, cex = gSize, col = gCol)


#######################
# Analysis of results #
#######################
# Basic information on input data can be accessed using the geneData function. The number of annotated
# genes, the number of significant genes (if it is the case), the minimal size of a GO category as well as the
# number of GO categories which have at least one signficant gene annotated are listed:
geneData(resultFisher)
# Annotated Significant    NodeSize    SigTerms
# 17944         677           5         253
# generate a summary of the enrichment analysis

results.table <- GenTable(GOdata.MF, resultFisher, topNodes = length(resultFisher@score))
# How many GO terms were tested?
dim(results.table)[1]
# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
results.table.bh <- results.table[which(p.adjust(results.table[,"result1"],method="BH")<=0.05),]
# How many terms are enriched?
dim(results.table.bh)[1]
# What are the top ten terms?
results.table.bh[1:10,]
# Get genes in most sig terms
# Get all the genes annotated to the top one GO term:
GOid.of.interest = results.table.bh[1,"GO.ID"]
all.term.genes = genesInTerm(GOdata.BP,GOid.of.interest)[[1]]
# Which of these genes is in the bicluster?
genes.of.interest <- intersect(genes,all.term.genes)

# generate a summary of the enrichment analysis
GenTable(GOdata, classicFisher = resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
#showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 50, useInfo = 'all')
printGraph(GOdata.MF, resultFisher, firstSigNodes = 20, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
if(is.na(sigNo<-table(score(result)<0.05)[2])){next}
showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
mtext(on, line=-1)

##############
# Automation #
##############
tgo.fisher<-function(genes)
{
    # genes of interest
    # head(genes)
    # limit genes to those already in university
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    # str(geneList)
    
    enrich<-list(MF=list(), BP=list(),CC=list())
    for(on in c("MF","BP","CC"))
    {
        print(on)
        # Make topGO object
        GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        # fisher test
        result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        results.table <- GenTable(GOdata, result, topNodes = length(result@score))
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
        results.table.bh <- results.table[which(p.adjust(results.table[,"result1"],method="BH")<=0.05),]
        enrich[[on]]<-list(GOproject=geneData(result),sig.bh=results.table.bh)
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-table(score(result)<0.05)[2])){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
    return(enrich)
}

