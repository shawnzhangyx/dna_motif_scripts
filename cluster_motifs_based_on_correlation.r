##################################
## This script will take a correlation matrix, and perform hierarchical clustering 
## and then cut the dendegram into clusters using a certain cutoff
##################################


inputfile = "~/software/meme/motifs/jaspar.core.vertebrate.motif.cor.csv"
mat = read.csv(inputfile,row.names=1) 
#mat.cut = t(mat[1:10,1:10])
colnames(mat) = rownames(mat)
#dis = as.dist(1-mat.cut)
dis = as.dist(1-t(mat))
hc = hclust(dis)
## use correlation 
classes = cutree(hc,h=0.2) 
class.sorted = sort(classes)
write.table(class.sorted, "~/software/meme/motifs/jaspar.core.vertebrate.motif_clusters.csv",quote=F,col.names=F,sep=',')


