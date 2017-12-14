glen<-unlist(geneset$setsize)

aa<-apply(tScores_mat,1,max)
bb<-apply(tScores_mat,1,min)

an1 = sqrt(2*log(glen))
top = log(4*pi)+log(log(glen))
bottom = 2*log(glen)
bn1 = an1 *(1-0.5*top/bottom)

newt<-aa
for ( i in 1: length(aa) ) {
  if ( abs(aa[i])< abs(bb[i]) ) { newt[i]<-bb[i] }
}


aa<-apply(tControl_mat,1,max)
bb<-apply(tControl_mat,1,min)


newc<-aa
for ( i in 1: length(aa) ) {
  if ( abs(aa[i])< abs(bb[i]) ) { newc[i]<-bb[i] }
}











mix.obj<-function(p,x,an, bn)
{
  z1<-(x-bn-p[2])*an*p[3]
  z2<-(x+bn+p[4])*an*p[5]
  e<-(p[1]*an*p[3])*exp(-z1-exp(-z1))+((1-p[1])*an*p[5])*exp(z2-exp(z2))
  if (any(e<=0)) Inf else -sum(log(e))
}
# mix.obj(p = p0, x = newc, an = an1, bn = bn1)

pp<-newc[newc>0]
nn<-newc[newc<0]

p0<-c(p=0.5,u1=1,s1=0.5 ,u2=1,s2=0.5)

lmix<-deriv(
  ~-log((p*an*s1)*exp(-((x-bn-u1)*an*s1)-exp(-(x-bn-u1)*an*s1))+((1-p)*an*s2)*exp(((x+bn+u2)*an*s2)-exp((x+bn+u2)*an*s2))),
  c("p","u1","s1","u2","s2"),
  function(x,an,bn,p,u1,s1,u2,s2) NULL)

mix.gr<-function(p,x,an,bn) {
  u1<-p[2]
  s1<-p[3]
  u2<-p[4]
  s2<-p[5]
  p<-p[1]
  colSums(attr(lmix(x,an,bn,p,u1,s1,u2,s2),"gradient"))}
# mix.gr(p = p0, x = newc, an = an1, bn = bn1)

aa<-optim(p0,mix.obj,mix.gr,x=newc,an=an1, bn=bn1, method="BFGS")

par<-aa$par
