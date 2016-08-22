
setwd("/Users/Doug/Downloads/SMC_vcf")


#!!!!!INDELS!!!!!

#LOAD VARIANT CONTROL FILE 1 
bp=read.table("SMC1.vcf")

#Extract relevent features
SMsub<- subset(bp, 
					V1=="ChLG3" 
					, select = c(V1,V2,V4,V5,V8,V9,V10))

###V1:Chromosome
###V2:Position on Chromosome
###V4:Reference Allele
###V5:Alternative Allele
###V8:Seqencing Info
###V9:penalized liklihood (not needed)
###V10:uncorrected sequencing info (also not needed) 


#Extracting relevent "Sequencing Info"
library(reshape2)

#we need the DP4 values: the Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles.

#Splitting V8 and adding dummy column names 
V8=colsplit(SMsub$V8,";",c("BlankM","SMC1M","SMC2M","ABM","GAM","NDGM","RajM","??","???"))
SMpracTotal=cbind(SMsub,V8)

#extracting only indel information for now.
INDELS=subset(SMpracTotal,BlankM=="INDEL")

#Splitting DP4 into 4 different columns
V8=colsplit(INDELS$RajM,",",c("?","??","???","????"))

#getting "forward ref alleles" as a numerical value without labels
V9=colsplit(V8[,1],"=",c("!","!!"))
yep=cbind(V9,V8)
yep[3]=NULL
yep[1]=NULL

#Combining all refferences alleles (a) and total alleles (b) to treat as a binomial response  
a=(yep[,1]+yep[,2])
b=(yep[,3]+yep[,4]+yep[,1]+yep[,2])
R=cbind(a,b)

#snps position information
position=INDELS$V2


#General Additive Binomial Fit of Ref Allele/Non-Ref Allele on chromosomal position  
library(mgcv)
x=bam(R~s(position),family=binomial)


#Plot of how the frequency of Ref to non-Ref Changes along the chromosome
plot(x, se = 1, seWithMean = TRUE, rug = TRUE, shift = mean(predict(x)),
        trans = function(x){exp(x)/(1+exp(x))},ylim=c(-10,2))
        
#Doesn't explain much        
        
        
px=predict(x)

#Saved for later
SM1=cbind(R,position)




###########################
#Now we are going to combine the high bulk (SMC1) and low bulk (SMC2) in one analysis
###########################


bp=read.table("SMC1.vcf")

SMsub<- subset(bp, 
					V1=="ChLG3" 
					, select = c(V1,V2,V4,V5,V8,V9,V10))

library(reshape2)

V8=colsplit(SMsub$V8,";",c("BlankM","SMC1M","SMC2M","ABM","GAM","NDGM","RajM","??","???"))
SMpracTotal=cbind(SMsub,V8)
INDELS=subset(SMpracTotal,BlankM=="INDEL")
V8=colsplit(INDELS$RajM,",",c("?","??","???","????"))
V9=colsplit(V8[,1],"=",c("!","!!"))

yep=cbind(V9,V8)

yep[3]=NULL
yep[1]=NULL

a=(yep[,1]+yep[,2])
b=(yep[,3]+yep[,4]+yep[,1]+yep[,2])
R=cbind(a,b)

position=INDELS$V2


#Same as the above example except we added the column "1" to distinguish these values from the 2nd bulk
SM1=cbind(R,position)
SM1=cbind(SM1,rep(1))




#Obtaining values for bulk 2
bp=read.table("SMC2.vcf")

SMsub<- subset(bp, 
					V1=="ChLG3" 
					, select = c(V1,V2,V4,V5,V8,V9,V10))

V8=colsplit(SMsub$V8,";",c("BlankM","SMC1M","SMC2M","ABM","GAM","NDGM","RajM","??","???"))
SMpracTotal=cbind(SMsub,V8)
INDELS=subset(SMpracTotal,BlankM=="INDEL")
V8=colsplit(INDELS$RajM,",",c("?","??","???","????"))
V9=colsplit(V8[,1],"=",c("!","!!"))

yep=cbind(V9,V8)

yep[3]=NULL
yep[1]=NULL

a=(yep[,1]+yep[,2])
b=(yep[,3]+yep[,4]+yep[,1]+yep[,2])
R=cbind(a,b)


position=INDELS$V2


#Values Saved and labeled 
SM2=cbind(R,position)
SM2=cbind(SM2,rep(2))




#combining the values from each bulk into one variable 
SMT=rbind(SM1,SM2)

#Separating Data into binomial response (R), chr location (position), and Bulk type  
R=cbind(SMT[,1],SMT[,2])
position=SMT[,3]
BULK=SMT[,4]



library(mgcv)

#modeling everything together regardless of bulk type

x21=bam(R~s(position,k=60), family=binomial)
plot(x21, se = 2, shade=TRUE, seWithMean = TRUE, rug = TRUE, shift = mean(predict(x21)),
        trans = function(x){exp(x)/(1+exp(x))},ylim=c(-10,5))


#because the high and low bulks are being modeled together the avg predicted allele freq across the chromosome should be .25

px21=predict(x21)
trans = function(x){exp(x)/(1+exp(x))}  #converting predicted values to the original scale
mean(trans(px21)) # mean should be 0.25, is 0.24485


#modeling bulks separately 

x2=bam(R~s(position,k=60)+s(position,by=BULK,k=60), family=binomial)
plot(x2, se = 2,shade=TRUE, seWithMean = TRUE, rug = TRUE, shift = mean(predict(x2)),
        trans = function(x){exp(x)/(1+exp(x))},ylim=c(-10,5))

#But what we really want is the maximum difference(s) between bulk 1 and bulk 2 
       
       
#####Differences
subsample=cbind(x2$model$position,x2$fitted.values,x2$model$BULK)
subsample1=subset(subsample,subsample[,3]==1)
subsample2=subset(subsample,subsample[,3]==2)
colnames(subsample1) <- c(paste("V",1:3))
colnames(subsample2) <- c(paste("V",1:3))
comsubsample=merge(subsample1,subsample2,by="V 1")

#Max difference between bulks 
max(comsubsample[,2]-comsubsample[,4])       

plot(spline(comsubsample[,1],(comsubsample[,2]-comsubsample[,4])))

#We now have the max difference between bulks. Now we need to build confidence intervals around that position

     

###Bootstrap
randomRows = function(df,n){
   return(df[sample(replace=TRUE,nrow(df),n),])
	}	
	
R=cbind(SMT[,1],SMT[,2])
position=SMT[,3]
BULK=SMT[,4]
	
BootX2=cbind(R,position,BULK)

nreps <- 5 											#Change me! 
nully=numeric(nreps+1)
nullx=numeric(nreps+1)

for (i in 2:nreps) {
	
	BootX21=randomRows(BootX2,length(BootX2[,1]))
	a=BootX21[,1]
	b=BootX21[,2]
	R=cbind(a,b)
	position=BootX21[,3]
	BULK=BootX21[,4]	
	x2r=bam(R~s(position,k=60)+s(position,by=BULK,k=60), family=binomial)
	
subsample=cbind(x2r$model$position,x2r$fitted.values,x2r$model$BULK)

subsample1=subset(subsample,subsample[,3]==1)
 
subsample2=subset(subsample,subsample[,3]==2)
        
colnames(subsample1) <- c(paste("V",1:3))

colnames(subsample2) <- c(paste("V",1:3))

comsubsample=merge(subsample1,subsample2,by="V 1")

comdiff=comsubsample[,2]-comsubsample[,4]

comsubsample=(cbind(comsubsample,comdiff))

bootmax=order(comsubsample$comdiff,decreasing=T)[1]

rowbootmax=comsubsample[bootmax,]

nully[i]=rowbootmax$comdiff
nullx[i]=rowbootmax[,1]

print(i)

	}

plot(spline(comsubsample[,1],(comsubsample[,2]-comsubsample[,4])/nreps))

print(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)

null1k=cbind(nullx,nully)

write.csv(null1k, file="/users/Doug/Desktop/null1k_2.csv")



#########SNPS!!!!!!!###########


bp=read.table("SMC1.vcf")

SMsub<- subset(bp, 
					V1=="ChLG3" 
					, select = c(V1,V2,V4,V5,V8,V9,V10))


library(reshape2)

V8=colsplit(SMsub$V8,";",c("BlankM","SMC1M","SMC2M","ABM","GAM","NDGM","RajM","??","???"))
SMpracTotal=cbind(SMsub,V8)
INDELS=subset(SMpracTotal,BlankM!="INDEL")


V7=colsplit(INDELS$SMC2M,"=",c("HK1","HK2"))
V7.7=cbind(INDELS,V7)

RPB=subset(V7.7,HK1=="RPB")
RPBnot=subset(V7.7,HK1!="RPB")


RPBsplit=colsplit(RPB$NDGM,",",c("?","??","???","????"))
RPBsplitFinal=colsplit(RPBsplit[,1],"=",c("!","!!"))
yep=cbind(RPBsplitFinal,RPBsplit)
yep[3]=NULL
yep[1]=NULL
a=(yep[,1]+yep[,2])
b=(yep[,3]+yep[,4]+yep[,1]+yep[,2])
RPBpos=RPB$V2
RPBcount=cbind(RPBpos,a,b)
RPBcount1=subset(RPBcount,b<=200)



RPBnotsplit=colsplit(RPBnot$GAM,",",c("?","??","???","????"))
RPBnotsplitFinal=colsplit(RPBnotsplit[,1],"=",c("!","!!"))
yep1=cbind(RPBnotsplitFinal,RPBnotsplit)
yep1[3]=NULL
yep1[1]=NULL
a1=(yep1[,1]+yep1[,2])
b1=(yep1[,3]+yep1[,4]+yep1[,1]+yep1[,2])
RPBnotpos=RPBnot$V2
RPBnotcount=cbind(RPBnotpos,a1,b1)
RPBnotcount1=subset(RPBnotcount,b1<=200)

SNPtotal=rbind(RPBcount1,RPBnotcount1)

R=cbind(SNPtotal[,2],SNPtotal[,3])

position=SNPtotal[,1]


SM1=cbind(R,position)
SM1=cbind(SM1,rep(1))





bp=read.table("SMC2.vcf")

SMsub<- subset(bp, 
					V1=="ChLG3" 
					, select = c(V1,V2,V4,V5,V8,V9,V10))


library(reshape2)

V8=colsplit(SMsub$V8,";",c("BlankM","SMC1M","SMC2M","ABM","GAM","NDGM","RajM","??","???"))
SMpracTotal=cbind(SMsub,V8)
INDELS=subset(SMpracTotal,BlankM!="INDEL")


V7=colsplit(INDELS$SMC2M,"=",c("HK1","HK2"))
V7.7=cbind(INDELS,V7)

RPB=subset(V7.7,HK1=="RPB")
RPBnot=subset(V7.7,HK1!="RPB")


RPBsplit=colsplit(RPB$NDGM,",",c("?","??","???","????"))
RPBsplitFinal=colsplit(RPBsplit[,1],"=",c("!","!!"))
yep=cbind(RPBsplitFinal,RPBsplit)
yep[3]=NULL
yep[1]=NULL
a=(yep[,1]+yep[,2])
b=(yep[,3]+yep[,4]+yep[,1]+yep[,2])
RPBpos=RPB$V2
RPBcount=cbind(RPBpos,a,b)
RPBcount1=subset(RPBcount,b<=200)



RPBnotsplit=colsplit(RPBnot$GAM,",",c("?","??","???","????"))
RPBnotsplitFinal=colsplit(RPBnotsplit[,1],"=",c("!","!!"))
yep1=cbind(RPBnotsplitFinal,RPBnotsplit)
yep1[3]=NULL
yep1[1]=NULL
a1=(yep1[,1]+yep1[,2])
b1=(yep1[,3]+yep1[,4]+yep1[,1]+yep1[,2])
RPBnotpos=RPBnot$V2
RPBnotcount=cbind(RPBnotpos,a1,b1)
RPBnotcount1=subset(RPBnotcount,b1<=200)

SNPtotal=rbind(RPBcount1,RPBnotcount1)

R=cbind(SNPtotal[,2],SNPtotal[,3])

position=SNPtotal[,1]



SM2=cbind(R,position)
SM2=cbind(SM2,rep(2))


SMT=rbind(SM1,SM2)
R=cbind(SMT[,1],SMT[,2])
position=SMT[,3]
BULK=SMT[,4]

x2=bam(R~s(position,k=60)+s(position,by=BULK,k=60), family=binomial)

plot(x2, se = 2,shade=TRUE, seWithMean = TRUE, rug = TRUE, shift = mean(predict(x2)),
        trans = function(x){exp(x)/(1+exp(x))},ylim=c(-10,5),)
        abline(h=0.25,col="red")
        abline(h=0.5,col="blue")
        abline(h=0.75,col="green")
        




head(SMT)

SMT95=(subset(bp,V2>38312931,V2< 38182988))







