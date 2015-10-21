library('glmnet')

x=read.table("glmnet_fit_inputs.txt",sep=",",col.names=c("x1","x2","x3","x4","x5","x6","x7","y"))
output<-as.matrix(data.frame(x$x1,x$x2,x$x3,x$x4,x$x5,x$x6,x$x7))
fit1<-glmnet(output, x$y, nlambda = 1, lambda = 20, alpha = 1)
print(coef(fit1))

cfs0<-coef(fit1)["(Intercept)",1]
cfs1<-coef(fit1)["x.x1",1]
cfs2<-coef(fit1)["x.x2",1]
cfs3<-coef(fit1)["x.x3",1]
cfs4<-coef(fit1)["x.x4",1]
cfs5<-coef(fit1)["x.x5",1]
cfs6<-coef(fit1)["x.x6",1]
cfs7<-coef(fit1)["x.x7",1]

cfs<-data.frame(cfs0,cfs1,cfs2,cfs3,cfs4,cfs5,cfs6,cfs7)
print(cfs)
write.table(cfs,"glmnet_fit_outputs.txt",sep=" ",row.names=FALSE,col.names=FALSE)
