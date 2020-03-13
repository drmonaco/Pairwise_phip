fread()

Pairwise  = function(full,tweak,bins){
  data = full
  data_dec = data[order(-data[,1]),, drop  = FALSE] # order the pairwsie points by the control (x) condition in decreasing order
  pred_1000 = rlm(data_dec[1:1000,2]~data_dec[1:1000,1]) # take the top 500 points (of the case condition) and get residual from that
  eq_after = full[,1]* coef(pred_1000)[2]+coef(pred_1000)[1]
  
  
  if(coef(pred_1000)[1]<0){
    pred_100 = rlm(data_dec[1:1000,2]~data_dec[1:1000,1]+0)
    eq_after = full[,1]* coef(pred_100)[1]
  }
  
  residual = full[,2]-eq_after # calculate the residual for all points in the pairwise comparison
  
  
  
  final = list()
  y=  bins
  sub_library = cbind(full,residual) # merge counts with residualss
  sub_library = sub_library[order(sub_library[,1]),,drop = FALSE]
  stdev  = c()
  med = c()
  
  for (R2 in 0:bins){
    Index = floor((1+(dim(sub_library)[1]/(y+1))*R2):(dim(sub_library)[1]/(y+1)+dim(sub_library)[1]/(y+1)*R2))
    Case_index = sub_library[(Index),]
    stdev[R2+1] = sd(Case_index[,3])
    med[R2+1] = median(Case_index[,1])
  }
  if(max(med[1:(bins-1)])!=0){
    eq = rlm(stdev[1:(bins-1)]~med[1:(bins-1)])
    x = (sub_library[,1]*coef(eq)[2]+coef(eq)[1])
    x[x<tweak] = tweak
    z_scores = sub_library[,3]/x
    
    
    if(coef(eq)[1]<0){
      eq = rlm(stdev[1:(bins-1)]~med[1:(bins-1)]+0)
      x = (sub_library[,3]*coef(eq)[2]+coef(eq)[1])
      x[x<tweak] = tweak
      z_scores = sub_library[,3]/x
    }
  }else {
    z_scores =  sub_library[,3]*0
  }
  final = cbind(sub_library,z_scores)
  return(final)
}

input = commandArgs(TRUE)


Utah = data.frame(fread(input[[1]]),row.names = 1)
colnames(Utah) = substring(colnames(Utah),2)

Utah2 = Utah
Utah = Utah2
colnames(Utah) = gsub("(20A20G.1).*","",colnames(Utah))
colnames(Utah) = paste0(colnames(Utah),"20A20G.1")

# colnames(Utah) = gsub("\\.","-",colnames(Utah))

input2 = data.frame(fread(input[[2]],header = FALSE))
input2 = input2[,c(2,1)]
tweak = 10
bins = 100
order = rownames(Utah)
out_path = input[[3]]


final_pairwise = data.frame(matrix(NA,nrow = dim(Utah)[1],ncol = dim(input2)[1]),row.names = row.names(Utah))
for(R in 1:dim(input2)[1]){
  
  ###
  pair = input2[R,]
  pair_name_1 = paste(c(pair[1],pair[2]),collapse = ":")
  pair = c(grep(pair[1],colnames(Utah)),  grep(pair[2],colnames(Utah)))
  full = Utah[,pair]
  if(sum(full[,1])<=1000 | sum(full[,2])<=1000){
    full = Utah[,pair]+1*rnorm(dim(full)[1])
  }
  final_up = Pairwise(full,tweak,bins)
  colnames(final_up) = c("before","after","residual","z_scores")
  final_up = final_up[match(order,rownames(final_up)),]
  
  
  ###
  pair = input2[R,]
  pair_name_2 = paste(c(pair[2],pair[1]),collapse = ":")
  pair = c(grep(pair[2],colnames(Utah)),  grep(pair[1],colnames(Utah)))
  full = Utah[,pair]
  if(sum(full[,1])<=1000 | sum(full[,2])<=1000){
    full = Utah[,pair]+1*rnorm(dim(full)[1])
  }
  final_down = Pairwise(full,tweak,bins)
  colnames(final_down) = c("before","after","residual","z_scores")
  final_down = final_down[match(order,rownames(final_down)),]
  
  
  final_pairwise[,R] = final_up$z_scores
  colnames(final_pairwise)[R] = pair_name_2
  final_pairwise[,R+dim(input2)[1]] = final_down$z_scores
  colnames(final_pairwise)[R+dim(input2)[1]] = pair_name_1
  print(pair_name_1)
  print(pair_name_2)
  full_plot = cbind(full,final_pairwise[,c(R,R+dim(input2)[1])],"black")
  full_plot[,6] = "black"
  full_plot[full_plot[,3]>=7,6] = "red"
  full_plot[full_plot[,4]>=7,6] = "blue"
  
  pair = input2[R,]
  if (input[[4]]=="yes"){
    full_matrix = cbind(full_plot,substr(rownames(full_plot),1,10))
    a = ggplot(full_matrix)+geom_count(aes(x = full_matrix[,1],y = full_matrix[,2],color = full_matrix[,6]))+facet_wrap(~full_matrix[,7],ncol = 2)+
      scale_x_log10()+scale_y_log10() + scale_colour_manual(values=c("black","red","blue"))+
      theme(legend.position = "none")+xlab(colnames(full_plot)[1])+ylab(colnames(full_plot)[2])
    b = ggplot(full_plot)+geom_count(aes(x = full_plot[,1],y = full_plot[,2],color = full_plot[,6]))+
      scale_x_log10()+scale_y_log10() + scale_colour_manual(values=c("black","red","blue"))+
      theme(legend.position = "none")+xlab(colnames(full_plot)[1])+ylab(colnames(full_plot)[2])
    ggsave(paste0(out_path,pair_name_1,".jpg"),arrangeGrob(b,a),height = 12,width = 7)
    
  } else{
    b = ggplot(full_plot)+geom_count(aes(x = full_plot[,1],y = full_plot[,2],color = full_plot[,6]))+
      scale_x_log10()+scale_y_log10()  + scale_colour_manual(values=c("black","red","blue"))+
      theme(legend.position = "none")+xlab(as.character(pair[1]))+ylab(as.character(pair[2]))
    ggsave(paste0(out_path,pair_name_1,".jpg"),b)
    
  }
  
}



rownames(final_pairwise) = order
fwrite(as.data.frame(round(final_pairwise,2)),paste0(out_path,"/pairwise_enrichments.csv"),row.names = TRUE)
