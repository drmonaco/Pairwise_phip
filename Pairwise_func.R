library(data.table)
library(tidyverse)
library(MASS)
library(OneR)

data = fread("pad_test_parwise.csv")
data_2 = data %>% select(ID,Unmod,PAD4)
#Pairwise  = function(full,tweak,bins){
data_2 = data_2 %>% arrange(-Unmod)
m1 <- rlm(data_2[1:500,3] ~ data_2[1:500,2], data = data_2[1:500,])
data_2 = data_2 %>% mutate(residual = data_2[,3]-(data_2[,2]*m1$coefficients[2]+m1$coefficients[1]))
sd_data = data_2 #%>% filter(Unmod >1) %>% filter(PAD4 >1 )

final = list()
y=  bins
sub_library = sd_data %>% select(-ID)
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
final$color = final$z_scores>10

ggplot(final,aes(x = Unmod+1,y = PAD4+1, color=  color))+geom_count()+scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")
