setwd("/Users/jhjoo/Desktop/var/final.codes/")
load(file="data/IRI.652759.prep.RData")


#### preprocessing predefined function  ####
IRIReshape=function(IRI){
  iri=list()
  category=unique(IRI$CATEGORY)
  for(i in 1:length(category)){
    iri[[as.character(category[i])]] = IRI[which(IRI$CATEGORY == category[i]), ]
  }
  return(iri)
}
#softmax=function(vector){
#  return(exp(vector)/sum(exp(vector)))
#}
SKUAgg=function(IRI.table, topk=10){
  #############################################
  # Top-K market-share SKU as cate represneted
  #############################################
  # UNITS: unot sales of SKU
  # SALES: dollars sales of SKU
  # PRICE: price of SKU
  # PZ.REDUCT: price reduction rate of SKU
  # PR: promotion decision of item
  # AD: advertising class of item
  # DI: display class of item.
  ############################################
  #IRI.table=IRI.652759.tables$diapers
  ############################################
  # Filter out the SKU out of k
  ############################################
  category=IRI.table$CATEGORY[1]
  skus = unique(IRI.table$COLUPC)
  if (length(skus) > topk){
    threshold=sort(unique(IRI.table$MS), decreasing=T)[topk+1]
    IRI.table=IRI.table[which(IRI.table$MS > threshold), ]
    skus=unique(IRI.table$COLUPC)
  }
  ############################################
  # Make up the lost week
  ############################################
  # Filter out the items without records more than half the year
  # And makeup the the records which > 6mo but < 12mo
  ############################################
  dd=mm=0
  for(i in 1:length(skus)){
    sku.sample=IRI.table[which(IRI.table$COLUPC==skus[i]), ]
    if (nrow(sku.sample)<(length(WEEKS.range)/2)-1){
      IRI.table=IRI.table[-(which(IRI.table$COLUPC==skus[i])), ]
      #cat('[DROP] SKU:', as.character(sku.sample$CATEGORY[1]),'/', skus[i], ' Times:', nrow(sku.sample), '\n')
      dd=dd+1
    } else if (nrow(sku.sample)<length(WEEKS.range)){
      makeup=makeupProcessing(sku.sample)
      IRI.table=rbind(IRI.table, makeup)
      #cat('[MAKE] SKU:', skus[i], ' Times:', nrow(sku.sample), '\n')
      mm=mm+1
    }
  }
  cat('Category:', as.character(category), '\n')
  cat('#EXIST:', length(skus), '#DROP:', dd, '#MAKE:', mm, '\n') 
  
  ############################################
  IRI.table$PZ.REDUCT=(IRI.table$PRICE - IRI.table$FIXED.PRICE) / IRI.table$FIXED.PRICE
  IRI.table$CATEGORY=IRI.table$FIXED.PRICE=IRI.table$MS=NULL
  if (nrow(IRI.table)>0){
    IRI.table$COLUPC <- paste(category, IRI.table$COLUPC, sep='')
  }
  # reshape
  data = reshape(IRI.table, direction='wide', idvar='WEEK', timevar='COLUPC')
  return(data)
}
CateAgg=function(IRI.table){
  ##############################################
  # Aggregate SKUs into category level variable
  ##############################################
  # UNITS: total unit sales of category
  # SALES: total dollars sales of category
  # PRICE: average price of SKU
  # PRICE.W: annual market-share-Weighted price of SKU
  # PZ.REDUCT: average price reduction ratio of each SKU
  # PR: promotion rate of category 
  # AD: advertising rate of category
  # DI: display rate of category
  # 
  # P.S. Assume SKUs are available in all year
  ############################################
  
  ############################################
  IRI.table$PRICE.W=IRI.table$PRICE * IRI.table$MS
  IRI.table$PZ.REDUCT= (IRI.table$PRICE - IRI.table$FIXED.PRICE) / IRI.table$FIXED.PRICE
  IRI.table$AD[which(IRI.table$AD > 0)] = 1
  IRI.table$DI[which(IRI.table$DI > 0)] = 1
  nsku = length(unique(IRI.table$COLUPC))
  category=IRI.table$CATEGORY[1]
  
  units = aggregate(IRI.table$UNITS~IRI.table$WEEK, IRI.table, sum)
  sales = aggregate(IRI.table$SALES~IRI.table$WEEK, IRI.table, sum)
  price = aggregate(IRI.table$PRICE~IRI.table$WEEK, IRI.table, mean)
  price.weighted = aggregate(IRI.table$PRICE.W~IRI.table$WEEK, IRI.table, sum)
  price.reduction = aggregate(IRI.table$PZ.REDUCT~IRI.table$WEEK, IRI.table, mean)
  promotion = aggregate(IRI.table$PR~IRI.table$WEEK, IRI.table, sum, na.rm=T)
  advertisement = aggregate(IRI.table$AD~IRI.table$WEEK, IRI.table, sum)
  display = aggregate(IRI.table$DI~IRI.table$WEEK, IRI.table, sum)
  
  #Merge
  data = Reduce(function(x, y) merge(x=x, y=y, by=c('IRI.table$WEEK')), 
                     list(units, sales, price, price.weighted, price.reduction, promotion, advertisement, display))
  colnames(data) <- c('WEEK', 'UNITS', 'SALES', 'PRICE', 'PRICE.W', 'PZ.REDUCT',  'PR.RATE', 'AD.RATE', 'DI.RATE')
  data$PR.RATE = data$PR.RATE / nsku
  data$AD.RATE = data$AD.RATE / nsku
  data$DI.RATE = data$DI.RATE / nsku
  colnames(data)[-1] <- paste(colnames(data)[-1], category, sep='.')
  
  return(data)
}
makeupProcessing=function(table){
  ############################################
  # Input the only records of this item.
  # Makeup the lost records with assumption.
  ############################################
  
  w.lost=setdiff(WEEKS.range, table$WEEK)
  makeup=do.call(rbind, replicate(length(w.lost), coredata(table[1, ]), simplify=FALSE))
  makeup$WEEK=w.lost
  makeup$UNITS=ifelse((mean(table$SALES) - 2*sd(table$SALES)<=0), 0, round(min(table$UNITS)))
  makeup$AD=min(Mode(table$AD))
  makeup$DI=min(Mode(table$DI))
  makeup$PR=min(Mode(table$PR))
  makeup$PRICE=ifelse(min(Mode(table$PR))==0, makeup$FIXED.PRICE, min(table$PRICE))
  makeup$SALES=makeup$UNITS * makeup$PRICE
  return(makeup)
}

#### category level ####
WEEKS.range=unique(IRI.652759.prep$WEEK)
IRI.652759.tables=IRIReshape(IRI.652759.prep)
IRI.652759.prep.cate=lapply(IRI.652759.tables, CateAgg)
save(IRI.652759.prep.cate, file="data/IRI.652759.prep.cate.RData")

#### sku level #####
WEEKS.range=unique(IRI.652759.prep$WEEK)
IRI.652759.tables=IRIReshape(IRI.652759.prep)
IRI.652759.prep.sku=lapply(IRI.652759.tables, SKUAgg, topk=10)
save(IRI.652759.prep.sku, file="data/IRI.652759.prep.sku.RData")

