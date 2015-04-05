library('glmnet')

x=read.table("glmnet_fit_inputs.txt",sep=",",col.names=c("x1","x2","x3","x4","x5","y"))
output<-as.matrix(data.frame(x$x1,x$x2,x$x3,x$x4,x$x5))
fit1<-glmnet(output,x$y)

cfs0<-coef(fit1)["(Intercept)",15]
cfs1<-coef(fit1)["x.x1",15]
cfs2<-coef(fit1)["x.x2",15]
cfs3<-coef(fit1)["x.x3",15]
cfs2<-coef(fit1)["x.x4",15]
cfs3<-coef(fit1)["x.x5",15]

cfs<-data.frame(cfs0,cfs1,cfs2,cfs3,cfs4,cfs5)
write.table(cfs,"glmnet_fit_outputs.txt",sep=" ",row.names=FALSE,col.names=FALSE)
