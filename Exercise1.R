# sampling t-distributed observations
returns=rt(10000,10)

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



#USING IDENTITY
# Computing the criterion function for each estimate...
output_1=c()
for (i in 1:length(para_list))
  {
    output_1=c(output_1,criterion_I(para_list[i],returns))
  }

# Finding the minimum value
minimum_1=para_list[which(output_1==min(output_1))]
print(minimum_1)


#... and then graph them!
pdf(file="t-returns_criterion_(W=I).pdf")
plot(para_list,output_1,type="l",xlab="nu",ylab="Criterion",main="W=I")
abline(v=minimum_1, col="purple", lty=2)
abline(v=10,col="red")
legend("topright", c("Estimated", "True Value"),col=c("purple","red"), lty=2:1, cex=1.0 )
dev.off()

# And charting the resulting density
pdf(file="t-returns_density_(W=I).pdf")
temp=density(scale(returns))
plot(temp, type="l", ylab="Density", main="t-returns density (W=I)")
lines(temp$x,dt(temp$x,minimum_1),col="red")
dev.off()

# USING COVARIANCE MATRIX
# Computing the criterion function for each estimate...
output_2=c()
for (i in 1:length(para_list))
{
  output_2=c(output_2,criterion_W(para_list[i],returns))
}

# Finding the minimum value
minimum_2=para_list[which(output_2==min(output_2))]
print(minimum_2)


#... and then graph them!
pdf(file="t-returns_criterion_(W=Sigma^-1).pdf")
plot(para_list,output_2,type="l",xlab="nu",
     ylab="Criterion",main="W=Cov(x)")
abline(v=minimum_2, col="purple",lty=2)
abline(v=10,col="red")
legend("topright", c("Estimated", "True Value"),col=c("purple","red"), lty=2:1, cex=1.0 )
dev.off()


# And charting the resulting density
pdf(file="t-returns_density_(W=Sigma^-1).pdf")
temp=density(scale(returns))
plot(temp, type="l", ylab="Density", main="t-returns density (W=Sigma^-1)")
lines(temp$x,dt(temp$x,minimum_2),col="red")
dev.off()
