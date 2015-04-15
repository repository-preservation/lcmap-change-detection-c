library('MASS')

x=read.table("robust_fit_inputs.txt",sep=",",col.names=c("x1","x2","x3","x4","y"))
output<-data.frame(x$y,x$x1,x$x2,x$x3,x$x4)
fit1<-rlm(x$y~x$x1+x$x2+x$x3+x$x4, data=output, psi=psi.bisquare, maxit = 5)

print(fit1)
cfs0=coef(summary(fit1))["(Intercept)","Value"]
cfs1=coef(summary(fit1))["x$x1","Value"]
cfs2=coef(summary(fit1))["x$x2","Value"]
cfs3=coef(summary(fit1))["x$x3","Value"]
cfs4=coef(summary(fit1))["x$x4","Value"]

cfs<-data.frame(cfs0,cfs1,cfs2,cfs3,cfs4)
write.table(cfs,"robust_fit_outputs.txt",sep=" ",row.names=FALSE,col.names=FALSE)
