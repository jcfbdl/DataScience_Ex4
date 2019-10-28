#########################################################################
# HEC Lausanne - MScF
# Data Science for Finance
# Exercise Session 2 - 27.09.2019
# Marceau Pierron, Taulant Ukshini, David Sasselli, Nora Koennyu
#########################################################################

# Setting Directory
setwd("/Users/marceau/Desktop/HEC_BA/MA1/DataScience_For_Finance/")

# Loading Data
y=read.delim("SP500_data.csv",sep=";",header=FALSE)

# Computing Log-Returns
y=y[,2]
y=diff(as.matrix(log(y)),1)
y=scale(y)

# Sampling t-distributed observations
#returns=rt(10000,10)
returns=y

# Creating a criterion function with W=I
# Step 1 Criterion
criterion_I<-function(para,x) {
  
  nu=para
  # 4th order moment condition
  cond1=mean(x^4)-((6)/(nu-4)+3)*(mean(x^2)^2)
  # 2nd order moment condition
  cond2=mean(x^2)-(nu/(nu-2))
  output_1=matrix(c(cond1,cond2),1,2)
  
  return(output_1%*%t(output_1))
  
}

# Creating a criterion function with the full W
# Step 2 Criterion
criterion_W<-function(para,x) {
  
  nu=para
  #4th order moment condition
  cond1=mean(x^4)-((6)/(nu-4)+3)*(mean(x^2)^2)
  #2nd order moment condition
  cond2=mean(x^2)-(nu/(nu-2))
  W=cov(cbind(x^2/sqrt(length(x)),x^4/sqrt(length(x))))
  output_2=matrix(c(cond1,cond2),1,2)
  
  return(output_2%*%solve(W)%*%t(output_2))
  
}

#########################################################################
# For W = I
#########################################################################

# Setting a potential list of candidate estimates
para_list_1=seq(4.1,30,by=0.01)
para_list_12=seq(4.1,7,by=0.01)
# Computing the criterion function for each estimate...
output_1=c()

for (i in 1:length(para_list_12)) {
  
  output_1=c(output_1,criterion_I(para_list_12[i],returns))
  
}

#... and then graph them:
plot(para_list_12,output_1,type="l",xlab="nu",ylab="Criterion",main="W=I")

# Finding the minimum value
minimum_1=para_list_1[which(output_1==min(output_1))]
print(minimum_1)

# And charting the resulting density
temp=density(returns)
plot(temp,type="l",main="")

lines(temp$x,dt(temp$x,5),col="red")

#########################################################################
# For W = W
#########################################################################

# Setting a potential list of candidate estimates
para_list_2=seq(4.1,30,by=0.01)
para_list_22=seq(4.1,7,by=0.01)
# Computing the criterion function for each estimate...
output_22=c()
output_2=c()

for (i in 1:length(para_list_22)) {
  
  output_22=c(output_22,criterion_W(para_list_22[i],returns))
  
}

#... and then graph them:
plot(para_list_22,output_22,type="l",xlab="nu",ylab="Criterion",main="W=W")

# Finding the minimum value
for (i in 1:length(para_list_2)) {
  
  output_2=c(output_2,criterion_W(para_list_2[i],returns))
  
}
minimum_2=para_list_2[which(output_2==min(output_2))]
print(minimum_2)

# And charting the resulting density
temp=density(scale(returns))
test_mean=mean(returns)
print(test_mean)
test_sd=sqrt(sum((returns-test_mean)^2)/(length(returns)-1))
print(test_sd)
plot(temp,type="l",main="")
lines(temp$x,dt(temp$x,4.2),col="blue")
lines(temp$x,dt(temp$x,5),col="red")
lines(temp$x,dt(temp$x,15),col="green")
lines(temp$x,dt(temp$x,100000),col="orange")

a=mean(returns^4)
b=mean(returns^2)
print(a)
print(b)

#eq = function(v){

#a=mean(x^4)
#b=mean(x^2)
#(a-(6/(v-4)+3)*b^2)^2+(b-v/(v-2))^2
#(mean(x^4)-((6)/(nu-4)+3)*(mean(x^2)^2))^2+(mean(x^2)-(nu/(nu-2)))^2
#(105.722215434505-(6/(v-4)+3)*2.08454826001654^2)^2+(2.08454826001654-v/(v-2))^2
#(40-(6/(v-4)+3)*2^2)^2+(2-v/(v-2))^2

#}

#curve(eq, from=0, to=16, xlab="v", ylab="crit")
