########################
## This sccript will take input from a MEME motif file
## and calculate the correlation between motif PWMs. 
## Timestamp: Yanxiao Zhang 2/20/2017 zhangyx.shawn@gmail.com
########################

input_filename = commandArgs(trailing=T)[1]


read_input_file_into_list = function(filename){
## read the motif PWMs into a list of data.frames
filename = "~/software/meme/motifs/vertebrate.motif"
lines = readLines(filename)
idx.start = which(startsWith(lines, "MOTIF"))

motifs = list()
for (idx in idx.start){
  incre = 2
  line = lines[idx+incre]
  name = strsplit(lines[idx]," ")[[1]][2]
  mat = as.numeric(unlist(strsplit(line,split='\t')))
  # if line is not empty, continue the loop
  while( line !=""){
  incre = incre + 1
  line = lines[idx+incre]
  prob = as.numeric(unlist(strsplit(line,split='\t')))
  mat = rbind(mat,prob)
  }
motifs[[name]] = mat
}
motifs
}

# motif1 = motifs[[1]]; motif2 = motifs[[2]]

compute_pwm_cor = function(motif1, motif2){
## calculate the average correlation between any two motif pairs. 
# make motif2 always the longer motif. 
if (length(motif1) > length(motif2)){
temp = motif2; motif2 = motif1;motif1 = temp
}
len1 = dim(motif1)[1]
len2 = dim(motif2)[1]
cor_list = NULL
# iterate through different location of the motif
for(iter in 1:len2){
motif2_comp = motif2[c(iter:len2,(1:iter-1)),]
# both sense and antisense.
  for (rev in c(FALSE,TRUE)){
    if (rev == TRUE){
      # reverse a motif. 
      motif2_comp = motif2_comp[len2:1,c(4,3,2,1)]
}
pho_list = NULL
for (y in 1:len1){
  pho = cor.test(motif1[y,],motif2_comp[y,])$estimate
  pho_list = c(pho_list,pho)
}
cor_list = c(cor_list,mean(pho_list))
}
}
cor= max(cor_list)
cor
}

compute_pwm_similarity_matrix = function(inputfilename,outmatrixname){
motifs = read_input_file_into_list(inputfilename)
len = length(motifs)
sim_mat = matrix(0,nrow=len,ncol=len)
colnames(sim_mat) = rownames(sim_mat) = names(motifs)
for (i in 1:(len-1)){
  for (j in (i+1):len){ 
  print(c(i,j))
    sim_mat[i,j] = compute_pwm_cor(motifs[[i]],motifs[[j]])
}
}
write.csv(sim_mat,outmatrixname,quote=F)
}

