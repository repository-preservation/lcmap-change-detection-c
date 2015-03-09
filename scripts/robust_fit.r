x=read.table("robust_fit.txt",sep=",",col.names=c("x1","x2","x3","x4","y"))
output<-data.frame(y,x1,x2,x3,x4)
fit1<-rlm(y~x1+x2+x3+x4, data=output,psi = psi.bisquare)
print(fit1)

cfs0=coef(summary(fit1))["(Intercept)","Value"]
cfs1=coef(summary(fit1))["x1","Value"]
cfs2=coef(summary(fit1))["x2","Value"]
cfs3=coef(summary(fit1))["x3","Value"]
cfs4=coef(summary(fit1))["x4","Value"]
print(cfs0)
print(cfs1)
print(cfs2)
print(cfs3)
print(cfs4)

cfs<-data.frame(cfs0,cfs1,cfs2,cfs3,cfs4)
print(cfs)
write.table(cfs,"robust_fit.txt",sep=" ",row.names=FALSE,col.names=FALSE)
