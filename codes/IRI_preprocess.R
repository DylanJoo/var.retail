setwd("/Users/jhjoo/Desktop/var/final.codes/")
load(file="data/IRI.652759.RData")

# UPC code convertion
upc.convert=function(lst.1, lst.2, lst.3, lst.4){
  upc=as.numeric(paste0(lst.1, lst.2, lst.3, lst.4))
  return(upc)
}
IRI.652759$COLUPC = upc.convert(IRI.652759$SY, IRI.652759$GE, IRI.652759$VEND, IRI.652759$ITEM)
IRI.652759$IRI_KEY = IRI.652759$SY = IRI.652759$GE = NULL
IRI.652759$VEND = IRI.652759$ITEM = IRI.652759$Year = IRI.652759$chain = NULL
colnames(IRI.652759)[which(colnames(IRI.652759) == 'category')] = 'CATEGORY'
colnames(IRI.652759)[which(colnames(IRI.652759) == 'DOLLARS')] = 'SALES'
colnames(IRI.652759)[which(colnames(IRI.652759) == 'F')] = 'AD'
colnames(IRI.652759)[which(colnames(IRI.652759) == 'D')] = 'DI'

# Get unit price
IRI.652759$PRICE = IRI.652759$SALES / IRI.652759$UNITS
IRI.652759[which(is.infinite(IRI.652759$PRICE)), ]$PRICE = IRI.652759[which(is.infinite(IRI.652759$PRICE)), ]$SALES

# Get fixed price
ModeorMean=function(x){
  # If no mode value 
  ret = ifelse(is.na(Mode(x)), mean(x), Mode(x))
  # If several mode value
  ret = ifelse(length(ret) > 1, max(ret), ret)
  return( ret )
}
temp=aggregate(IRI.652759$PRICE~IRI.652759$PR+IRI.652759$COLUPC, IRI.652759, ModeorMean)
category.fixed.price=aggregate(temp$`IRI.652759$PRICE`~temp$`IRI.652759$COLUPC`, temp, max)
IRI.652759 = merge(IRI.652759, category.fixed.price, by.x= c('COLUPC'), by.y=c('temp$`IRI.652759$COLUPC`'))
colnames(IRI.652759)[ncol(IRI.652759)] <- 'FIXED.PRICE'

# Get market-share of item by (e.g. BEER_1 annual sales / BEER annual total sales)
category.annual.SALES = aggregate(IRI.652759$SALES~IRI.652759$CATEGORY, IRI.652759, sum)
IRI.652759 = merge(IRI.652759, category.annual.SALES, by.x= c('CATEGORY'), by.y=c('IRI.652759$CATEGORY'))
colnames(IRI.652759)[ncol(IRI.652759)] <- 'SALES.CATE'

sku.annual.SALES = aggregate(IRI.652759$SALES~IRI.652759$COLUPC, IRI.652759, sum)
IRI.652759 = merge(IRI.652759, sku.annual.SALES, by.x= c('COLUPC'), by.y=c('IRI.652759$COLUPC'))
colnames(IRI.652759)[ncol(IRI.652759)] <- 'SALES.SKU'

IRI.652759$MS = IRI.652759$SALES.SKU / IRI.652759$SALES.CATE
IRI.652759$SALES.CATE = IRI.652759$SALES.SKU = NULL

# Encoding the categorical data format. e.g. Advertising(AD), Display(DI), 
IRI.652759$AD = as.numeric(as.character(factor(IRI.652759$AD, levels = c("A", "A+", "B", "C", "NONE"), labels = c(3, 4, 2, 1, 0))))
IRI.652759$DI = as.numeric(as.character(factor(IRI.652759$DI, levels = c(2, 1, 0), labels = c(2, 1, 0))))

# rename 
IRI.652759.prep=IRI.652759
save(IRI.652759.prep, file="data/IRI.652759.prep.RData")