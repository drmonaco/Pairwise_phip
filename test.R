
x = sd_data %>% mutate(test = bin(sd_data$Unmod,nbins = 50,method = "content")) %>% group_by(test) %>% summarise(sd(Unmod)) %>% 
  filter(`sd(Unmod)` !=0)

temp  = sd_data %>% mutate(x = ntile(temp$Unmod,10)) #%>% group_by(x) %>% summarise(sd(Unmod))# %>% 
filter(`sd(Unmod)` !=0)


  
  ntile_ <- function(x, n) {
    b <- x[!is.na(x)]
    q <- floor((n * (rank(b, ties.method = "first") - 1)/length(b)) + 1)
    d <- rep(NA, length(x))
    d[!is.na(x)] <- q
    return(d)
  }
  
  