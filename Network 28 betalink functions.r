##This file contains functions written by Poisot and colleagues to implement their analysis of dissimilarity of species interaction networks (Poisot et al 2012 Ecology Letters).
##Rabie's adjustments are noted in the code

library(vegan)

B01 = function(pm) with(pm,{2*(a+b+c)/(2*a+b+c) - 1})
#B02 = function(pm) B01(pm)
#B03 = function(pm) with(pm,{(b+c)/2})
#B04 = function(pm) with(pm,{(b+c)})
#B05 = function(pm) with(pm,{(((a+b+c)^2)/((a+b+c)^2-(2*b*c)))-1})
#B06 = function(pm) with(pm,{log(2*a+b+c)-((2*a*log(2))/(2*a+b+c))-(((a+b)*log(a+b)+(a+c)*log(a+c))/(2*a+b+c))})
#B07 = function(pm) exp(B06(pm))-1
#B08 = function(pm) with(pm,{(b+c)/(2*a+b+c)})
#B09 = function(pm) B08(pm)
#B10 = function(pm) with(pm,{a/(a+b+c)})
#B11 = function(pm) with(pm,{(2*a)/(2*a+b+c)})
#B12 = function(pm) with(pm,{(2*a+b+c)*(1-(a/(a+b+c)))})
#B13 = function(pm) with(pm,{min(b,c)/(max(b,c)+a)})
#B14 = function(pm) with(pm,{1-((a*(2*a+b+c))/(2*(a+b)*(a+c)))})
#B15 = function(pm) with(pm,{(b+c)/(a+b+c)})
#B16 = function(pm) B15(pm)
#B17 = function(pm) with(pm,{min(b,c)/(a+b+c)})
#B18 = function(pm) with(pm,{(b+c)/2})
#B19 = function(pm) with(pm,{(b*c+1)/(((a+b+c)^2-(a+b+c))/2)})
#B20 = function(pm) with(pm,{1-(2*a)/(2*a+b+c)/2})
#B21 = function(pm) with(pm,{a/(a+c)})
#B22 = function(pm) with(pm,{min(b,c)/(min(b,c)+a)})
#B23 = function(pm) with(pm,{(2*aba(b-c))/(2*a+b+c)})
#B24 = function(pm) with(pm,{1-(log((2*a+b+c)/(a+b+c))/log(2))})

##Rabie modified to make this function amenable to quantitative distance metrics.  
  ##metaweb is now sum of observations over all realizations, 
  ##co.oc is still the number of realizations for which each species pair co-occurs, and 
  ##null.template is now the expected number of observed interactions per co-occurence.
aggregate.metaweb = function(W,binary=F){
	Lo = unique(unlist(lapply(W,colnames)))
	Up = unique(unlist(lapply(W,rownames)))
	meta = matrix(0,ncol=length(Lo),nrow=length(Up))
	colnames(meta) = Lo
	rownames(meta) = Up
	co.oc = meta
	for(w in W){
  if(binary==T)	w[w>0] = 1   ##binarization is now optional.
		meta[rownames(w),colnames(w)] = meta[rownames(w),colnames(w)] + w
		co.oc[rownames(w),colnames(w)] = co.oc[rownames(w),colnames(w)] + 1
	}
	null.template = meta/co.oc
	null.template[is.nan(null.template)] = 0
#	meta[meta>0] = 1        ##remove binarization.
	return(list(web=meta,template=null.template, cooc = co.oc))
}

#now amenable to quantitative metrics
beta.os_prime = function(W,bf="jaccard"){
	metaweb = aggregate.metaweb(W)$web
	os_prime = NULL
	for(w in W){
		partitions = betalink(as.table(w),as.table(metaweb[rownames(w),colnames(w)]),bf)
		os_prime = c(os_prime,partitions$OS)
	}
	return(os_prime)
}


##Rabie reduced text size for row and column names.  I can imagine adjusting alpha values based on interaction freqencies to make this a quantitative visualization aid, but that's low priority for me.
betalink.plot = function(m1,m2,by.unique=TRUE){
	up.1 = rownames(m1)
	up.2 = rownames(m2)
	lo.1 = colnames(m1)
	lo.2 = colnames(m2)
	up.shared = up.1[up.1%in%up.2]
	lo.shared = lo.1[lo.1%in%lo.2]
	up.un.1 = up.1[!(up.1%in%up.shared)]
	up.un.2 = up.2[!(up.2%in%up.shared)]
	lo.un.1 = lo.1[!(lo.1%in%lo.shared)]
	lo.un.2 = lo.2[!(lo.2%in%lo.shared)]
	m = aggregate.metaweb(list(m1,m2))$web
	if(by.unique){
		oCol = rep(0,ncol(m))
		oCol[colnames(m)%in%lo.un.1] = c(1:length(lo.un.1))
		oCol[colnames(m)%in%lo.shared] = max(oCol)+c(1:length(lo.shared))
		oCol[colnames(m)%in%lo.un.2] = max(oCol)+c(1:length(lo.un.2))
		oRow = rep(0,nrow(m))
		oRow[rownames(m)%in%up.un.1] = c(1:length(up.un.1))
		oRow[rownames(m)%in%up.shared] = max(oRow)+c(1:length(up.shared))
		oRow[rownames(m)%in%up.un.2] = max(oRow)+c(1:length(up.un.2))
	} else {
		oCol = ncol(m)-rank(colSums(m),ties.method='r')+1
		oRow = nrow(m)-rank(rowSums(m),ties.method='r')+1
	}
	plot(0, pch=NA, xlim=c(-3,nrow(m)), ylim=c(-3,ncol(m)), asp=1, ylab='',xlab='', yaxt='n', xaxt='n', bty='n')
	text(y=-0.4,x=oRow+0.5,rownames(m),pos=2,srt=90,cex=0.5)
	text(x=-0.4,y=oCol-0.5,colnames(m),pos=2,cex=0.5)
	for(ro in 1:nrow(m)){
		current.up = rownames(m)[ro]
		for(co in 1:ncol(m)){
			current.low = colnames(m)[co]
			c.BG = 'lightgrey'
			if((current.up%in%up.un.1) | (current.low%in%lo.un.1)){c.BG = 'palegreen'}
			if((current.up%in%up.un.2) | (current.low%in%lo.un.2)){c.BG = 'skyblue'}
			if((current.up%in%up.shared) & (current.low%in%lo.shared)){
				if(m1[current.up,current.low]!=m2[current.up,current.low]){
					c.BG = 'orange'
				}
			}
			if(m[current.up,current.low] == 0){c.BG = 'white'}
			rect(oRow[ro]-0.9,oCol[co]-0.9,oRow[ro]-0.1,oCol[co]-0.1,col=c.BG,border=NA)
		}
	}
}

##Modified so that betalink.dist returns distance matrices for upper and lower trophic levels in the bipartite matrix.  
##added an option to accommodate quantitative metrics
##Also added an option for rectangular output of distance matrices
betalink.dist = function(W,triangular=F,bf="jaccard"){
	dWN = matrix(NA,ncol=length(W),nrow=length(W))
	colnames(dWN) = names(W)
	rownames(dWN) = names(W)
	diag(dWN)=0
	dOS = dWN
	dS_all.taxa = dWN
	dS_insects = dWN
	dS_plants = dWN
	dST = dWN
	dContrib = dWN
	for(i in 1:(length(W)-1)){  ##removed unnecessary c()
		for(j in (i+1):length(W)){  ##removed unnecessary c()
      	index=cbind(c(j,i),c(i,j))
#      	cat(names(W)[i],names(W)[j],"\n")
 			partition = betalink(W[[i]],W[[j]],bf)
			dWN[index]          = partition$WN
			dOS[index] 	        = partition$OS
			dS_all.taxa[index]	= partition$S
      dS_insects[index]   = partition$U
    	dS_plants[index]    = partition$L
			dST[index]          = partition$ST
			dContrib[index]     = partition$contrib
		}
	}
	distances = list(WN=dWN, OS=dOS, S_all.taxa=dS_all.taxa, S_insects=dS_insects ,S_plants=dS_plants, ST=dST, contrib=dContrib)
  if(triangular==T) distances = lapply(distances,as.dist) 
	return(distances)
}


##betalink is now a dispatch function that calls its binary version if the bf argument is itself a function, and calls its quantitative version if the bf argument is a character string
betalink=function(w1,w2,bf="jaccard"){
  if(is.character(bf)) betalink.q(w1,w2,bf=bf) else betalink.b(w1,w2,bf=bf)
}

#This used to be called betalink.  It calculates binary distance metrics.  If given a quantitative network it will convert it to binary with a warning.
#Rabie changed output of    		
  #beta_WN 
  #beta_ST 
	#b_contrib 
# for cases where there is no overlap between species composition of the two networks.  My reasoning is that the networks must be completely different so the dis-similarity values for these components must be 1.
betalink.b = function(w1,w2,bf=B01){
  if(any(!unique(c(w1,w2)) %in% c(0,1))) {
    warning("Quantitative matrix converted to binary matrix")
    w1[w1>0]=1
    w2[w2>0]=1 
    }
	pmb = function(A,B) list(b=sum(!(A %in% B)),c=sum(!(B %in% A)),a=sum(B %in% A))
	sp1 = list(top=rownames(w1),bottom=colnames(w1),all=unique(c(colnames(w1),rownames(w1))))
	sp2 = list(top=rownames(w2),bottom=colnames(w2),all=unique(c(colnames(w2),rownames(w2))))
	beta_U = bf(pmb(sp1$top,sp2$top))
	beta_L = bf(pmb(sp1$bottom,sp2$bottom))
	beta_S = bf(pmb(sp1$all,sp2$all))
	# Common species
	Csp = sp1$all[sp1$all %in% sp2$all]
	CUsp = sp1$top[sp1$top %in% sp2$top]
	CLsp = sp1$bottom[sp1$bottom %in% sp2$bottom]
	if((length(CUsp)>0) & (length(CLsp)>0))
	{
		w1Con = w1[CUsp,CLsp]
		w2Con = w2[CUsp,CLsp]
		nCon = sum((w1Con == w2Con) & (w1Con == 1))
		pmBos = list(b=sum(w1Con)-nCon,c=sum(w2Con)-nCon,a=nCon)
		pmBwn = list(b=sum(w1)-nCon,c=sum(w2)-nCon,a=nCon)
		beta_OS = bf(pmBos)
		beta_WN = bf(pmBwn)
		if(is.na(beta_OS)) beta_OS = 0
		if(is.na(beta_WN)) beta_WN = 0
		beta_ST = beta_WN - beta_OS
		if(beta_WN > 0){
			b_contrib = beta_ST / beta_WN
		} else {
			b_contrib = 0
		}
	} else {
		beta_WN = 1
		beta_OS = 0
		beta_ST = 1
		b_contrib = 1
	}
	return(list(U = beta_U, L = beta_L, S = beta_S, OS = beta_OS, WN = beta_WN, ST = beta_ST, contrib = b_contrib))
}

#if one data frame is a subset of a second data frame, the merge function will return the larger data frame.  This function forces the merge function to preserve all rows from both data sets
vec2data.frame=function(x,y,all=T){
  temp=merge(t(x),t(y),all=all)
  if(nrow(temp)<sum(min(nrow(x),1),min(nrow(y),1))){
  x=c(a = -1,x)
  y=c(a = -2,y)
  temp=merge(t(x),t(y),all=all)
  n=names(temp)
  temp=data.frame(temp[,-1])
  names(temp)=n[-1]}
  return(temp)}

#Quantitative version of the betalink function.    
#Treatment of beta_WN, beta_ST and contrib is the same as in betalink.b  See notes associated with that function.  
betalink.q = function(w1,w2,bf="jaccard"){
	sp1 = list(top=rowSums(w1),bottom=colSums(w1),all=c(rowSums(w1),colSums(w1)[!names(colSums(w1))%in%names(rowSums(w1))]))
	sp2 = list(top=rowSums(w2),bottom=colSums(w2),all=c(rowSums(w2),colSums(w2)[!names(colSums(w2))%in%names(rowSums(w2))]))
	topCom=vec2data.frame(sp1$top,sp2$top,all=T)
  bottomCom=vec2data.frame(sp1$bottom,sp2$bottom,all=T)
  allCom=vec2data.frame(sp1$all,sp2$all,all=T)  
 	topCom[is.na(topCom)]=0
  bottomCom[is.na(bottomCom)]=0
  allCom[is.na(allCom)]=0
	beta_U = vegdist(topCom,method=bf,diag=T,upper=T)
 	beta_L = vegdist(bottomCom,method=bf,diag=T,upper=T)
	beta_S = vegdist(allCom,method=bf,diag=T,upper=T)
	# Common species
	Csp = names(sp1$all)[names(sp1$all) %in% names(sp2$all)]
	CUsp = names(sp1$top)[names(sp1$top) %in% names(sp2$top)]
	CLsp =names(sp1$bottom)[names(sp1$bottom) %in% names(sp2$bottom)]
	if((length(CUsp)>0) & (length(CLsp)>0))
	{
		w1ConInts = as.data.frame(as.table(w1[CUsp,CLsp]))$Freq
		names(w1ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
		w2ConInts= as.data.frame(as.table(w2[CUsp,CLsp]))$Freq
		names(w2ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
    commonTaxaWeb=vec2data.frame(w1ConInts,w2ConInts,all=T)   
    w1=as.data.frame(w1)
    w1Ints = w1$Freq
    names(w1Ints) = paste(w1$Var1,w1$Var2,sep="-")
    w2=as.data.frame(w2)
    w2Ints = w2$Freq
    names(w2Ints) = paste(w2$Var1,w2$Var2,sep="-")
    ints=vec2data.frame(w1Ints,w2Ints,all=T)
    ints[is.na(ints)]=0  
    beta_WN = vegdist(ints,method=bf,diag=T,upper=T)
    if(all(rowSums(commonTaxaWeb)>0)){
      beta_OS = vegdist(commonTaxaWeb,method=bf,diag=T,upper=T)} else   {
        beta_OS = 1
        }
		beta_ST = beta_WN - beta_OS
		if(beta_WN > 0){
			b_contrib = beta_ST / beta_WN
		} else {
			b_contrib = 0
		}
	} else {
		beta_WN = 1
		beta_OS = 0
		beta_ST = 1
		b_contrib = 1 
	}
	return(list(U = beta_U[1], L = beta_L[1], S = beta_S[1], OS = beta_OS[1], WN = beta_WN[1], ST = beta_ST[1], contrib = b_contrib[1]))
}


##########################Null model analysis.  Rabie did not touch these.
#
#extract.localwebs = function(template,u,l,i){
#	ZDS = TRUE
#	while(ZDS){
#		w1 = template[sample(c(1:nrow(template)),u),sample(c(1:ncol(template)),l)]
#		w2 = template[sample(c(1:nrow(template)),u),sample(c(1:ncol(template)),l)]
#		## Link inconsistency
#		for(uu in 1:u) for(ll in 1:l){
#			if(w1[uu,ll] == 1){
#				if(rbinom(1,1,i) == 0){
#					if(runif(1,0,1)<0.5){
#						w1[uu,ll] = 0
#					} else {
#						w2[uu,ll] = 0
#					}
#				}
#			}
#		}
#		ZDS = sum(colSums(w1) == 0)+sum(rowSums(w1) == 0)+sum(colSums(w2) == 0)+sum(rowSums(w2) == 0) > 0
#	}
#	return(list(m1=w1,m2=w2))
#}
#
#
#null.from_template = function(w,template,reps){
#	temp = template[rownames(w),colnames(w)]
#	temp = temp / max(temp)
#	nulls_mw = null.metaweb(temp,reps)
#	nulls_de = null.degrees(temp,reps)
#	nulls_co = null.connectance(temp,reps)
#	return(I = nulls_co, II = nulls_de, MW = nulls_mw)
#}
#
#null.metaweb = function(template,reps)
#{
#	LNulls = NULL
#	n = 0
#	while(n < reps)
#	{
#		Void = TRUE
#		while(Void)
#		{
#			tNull = matrix(rbinom(prod(dim(template)),1,template),ncol=ncol(template),byrow=TRUE)
#			if ((sum(colSums(tNull) == 0)+sum(rowSums(tNull) == 0)) == 0) Void = FALSE
#		}
#		n = n+1
#		colnames(tNull) = colnames(template)
#		rownames(tNull) = rownames(template)
#		LNulls[[n]] = tNull
#	}
#	return(LNulls)
#}
#
#null.degrees = function(template,reps)
#{
#	template[template>0] <- 1
#	## Get proba
#	margin <- ifelse(ncol(template)<nrow(template),2,1)
#	currentweb <- matrix(0,ncol=ncol(template),nrow=nrow(template))
#	pc <- colMeans(template)
#	pr <- rowMeans(template)
#	if(margin==2)
#	{
#		for(i in 1:ncol(template))
#		{
#			currentweb[,i] <- (pc[i]+pr)/2
#		}
#	} else {
#		for(i in 1:nrow(template))
#		{
#			currentweb[i,] <- (pr[i]+pc)/2
#		}
#	}
#	#
#	LNulls = NULL
#	n = 0
#	while(n < reps)
#	{
#		Void = TRUE
#		while(Void)
#		{
#			tNull = apply(currentweb,margin,function(x)rbinom(length(x),1,x))
#		if ((sum(colSums(tNull) == 0)+sum(rowSums(tNull) == 0)) == 0) Void = FALSE
#		}
#		n = n+1
#		colnames(tNull) = colnames(template)
#		rownames(tNull) = rownames(template)
#		LNulls[[n]] = tNull
#	}
#	return(LNulls)	
#}
#
#null.connectance = function(template,reps)
#{
#	template[template>0] <- 1
#	LNulls = NULL
#	n = 0
#	while(n < reps)
#	{
#		Void = TRUE
#		while(Void)
#		{
#			tNull = matrix(sample(template),nrow=nrow(template))
#		if ((sum(colSums(tNull) == 0)+sum(rowSums(tNull) == 0)) == 0) Void = FALSE
#		}
#		n = n+1
#		colnames(tNull) = colnames(template)
#		rownames(tNull) = rownames(template)
#		LNulls[[n]] = tNull
#	}
#	return(LNulls)	
#}
#
#
#
#
#generate.metaweb = function(U,L,k){
#	Nzeroes = U*L - k
#	ZDS = TRUE
#	while(ZDS == TRUE){
#		Adj = sample(c(rep(0,Nzeroes),rep(1,k)),replace=FALSE)
#		W = matrix(Adj,nrow=U)
#		ZDS = sum(colSums(W) == 0)+sum(rowSums(W) == 0) > 0
#	}
#	colnames(W) = paste('l',c(1:ncol(W)), sep='')
#	rownames(W) = paste('u',c(1:nrow(W)), sep='')
#	return(W)
#}
#