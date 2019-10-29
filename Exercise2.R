x=read.delim("sp.csv",sep=";",header=FALSE)
x=diff(as.matrix(log(x)),1)

# sampling t-distributed observations
returns=(x-mean(x))/sd(x)
       
# Creating a criterion function with W=I
criterion_I<-function(para,x)
 {
   nu=para
   cond1=mean(x^4)-((6)/(nu-4)+3)*mean(x^2)^2
   cond2=mean(x^2)-(nu/(nu-2))
   # output=y^Ty
   output=matrix(c(cond1,cond2),1,2)
   return(output%*%t(output))
 }
       
# Creating a criterion function with the full W
criterion_W<-function(para,x)
 {
   nu=para
   cond1=mean(x^4)-((6)/(nu-4)+3)*mean(x^2)^2
   cond2=mean(x^2)-(nu/(nu-2))
   W=cov(cbind(x^4,x^2))
   output=matrix(c(cond1,cond2),1,2)
   return(output%*%solve(W)%*%t(output))
 }
       
# Setting a potential list of candidate estimates
para_list=seq(5,30,by=1)
     
       
       
#USING IDENTITY MATRIX
# Computing the criterion function for each estimate...
output_3=c()
for (i in 1:length(para_list))
 {
   output_3=c(output_3,criterion_I(para_list[i],returns))
 }
       
# Finding the minimum value
minimum_3=para_list[which(output_3==min(output_3))]
print(minimum_3)

#... and then graph them!
pdf(file="S&P500_returns_criterion_(W=I).pdf")
plot(para_list,output_3,type="l",xlab="nu",ylab="Criterion",main="W=I")
abline(v=minimum_3, col="purple",lty=2)
legend("bottomright", c("Estimated"),col=c("purple"), lty=2, cex=1.0 )
dev.off()

       
# And charting the resulting density
pdf(file="S&P500_returns_density_(W=I).pdf")
temp=density(scale(returns))
plot(temp, type="l", ylab="Density", main="S&P500 returns density (W=I)")
lines(temp$x,dt(temp$x,minimum_3),col="red")       
dev.off()


# USING COVARIANCE MATRIX
# Computing the criterion function for each estimate...
output_4=c()
for (i in 1:length(para_list))
 {
   output_4=c(output_4,criterion_W(para_list[i],returns))
 }
       
# Finding the minimum value
minimum_4=para_list[which(output_4==min(output_4))]
print(minimum_4)


#... and then graph them!
pdf(file="S&P500_returns_criterion_(W=Sigma^-1).pdf")
plot(para_list,output_4,type="l",xlab="nu",ylab="Criterion",main="W=Cov(x)")
abline(v=minimum_4, col="purple",lty=2)
legend("top", c("Estimated"),col=c("purple"), lty=2, cex=1.0 )
dev.off()


# And charting the resulting density
pdf(file="S&P500_returns_density_(W=Sigma^-1).pdf")
temp=density(scale(returns))
plot(temp, type="l", ylab="Density", main="S&P500 returns density (W=Sigma^-1)")
lines(temp$x,dt(temp$x,minimum_4),col="red")       
dev.off()