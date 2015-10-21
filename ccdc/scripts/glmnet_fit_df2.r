library('glmnet')

x=read.table("glmnet_fit_inputs.txt",sep=",",col.names=c("x1","y"))
output<-as.matrix(data.frame(x$x1))
fit1<-glmnet(output, x$y, nlambda = 1, lambda = 20, alpha = 1)
print(coef(fit1))

cfs0<-coef(fit1)["(Intercept)",1]
cfs1<-coef(fit1)["x.x1",1]

cfs<-data.frame(cfs0,cfs1)
print(cfs)
write.table(cfs,"glmnet_fit_outputs.txt",sep=" ",row.names=FALSE,col.names=FALSE)
