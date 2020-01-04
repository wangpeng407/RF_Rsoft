library(optparse)

option_list <- list(
  make_option(c("-i","--infile"), action="store", type="character", default=NULL, help="Input the data, otus, genus, or other functional abundance matirx."),
  make_option(c("-o", "--outdir"), action="store", default="./", type="character", help="The output dirctory, default is ./"),
  make_option(c("-t", "--train_set_proportion"), action="store", default=0.7, type="double", help="Set proportion of tranning dataset for models, default is 0.7"),
  make_option(c("-k", "--k_fold"), action="store", default=5, type="integer", help="K-fold cross validation, default is 5"),
  make_option(c("-m", "--mtry"), action="store", default=NULL, type="integer", help="mtry for random forest, if not set, cv will deterimine the best"),
  make_option(c("-n", "--ntree"), action="store", default=NULL, type="integer", help="ntree for random forest, if not set, cv will deterimine the best"),
  make_option(c("-p", "--top"), action="store", default=30, type="integer", help="Choose the top variable for vasulization of MeanDecreaseAccuracy and MeanDecreaseGini, default is 30")
)

opt <- parse_args(OptionParser(usage="%prog [options] file\n", option_list=option_list))

package_list <- c("reshape2","ggplot2","randomForest", "pROC")
for(pk in package_list){
  if(!suppressWarnings(suppressMessages(require(pk, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    stop("WARNING: Please install ", pk, "!\n")
  }else{
    cat("YES:", pk, "was succesfully installed and loaded!\n")
  }
}
cat('####################START:', date(),'####################\n')
###FUNCTIONS###
boxdata_get <- function(dt){
  gn <- names(dt)[-1]
  gs <- as.vector(dt$Class)
  temp.gs <- c()
  temp.value <- c()
  for(i in 1:length(gn)){
    idx <- which(gn[i] == gs)
    temp.gs <- c(temp.gs, gs[idx])
    temp.value <- c(temp.value, dt[idx, i+1])
  }
  temp.gs <- factor(temp.gs, levels = unique(dt$Class))
  bdt <- data.frame(Class = temp.gs, Prob = temp.value)
  return(bdt)
}

multiple_roc_plot <- function(responce, predictor, outdir = './'){
  cmbind <- as.data.frame(combn(as.vector(unique(responce)), 2))
  if(is.na(file.info(outdir)$isdir)) dir.create(outdir, recursive = TRUE)
  for(i in 1:ncol(cmbind)){
    cg <- cmbind[,i]
    outpdf <- paste0(outdir, '/', paste(cg, collapse = '-'), '.roc.pdf')
    pdf(outpdf, width = 8, height = 8)
    multiclass.roc(responce, predictor, percent = F, plot = T,
                   print.auc = T, auc.polygon=T, grid=c(0.1,0.1), 
                   legacy.axes = T, levels = cg,
                   grid.col=c("gray", "gray"), max.auc.polygon=T,
                   auc.polygon.col="lightgreen")
    dev.off()
  }
}

if(is.na(file.info(opt$outdir)$isdir)) dir.create(opt$outdir, recursive = TRUE)
if(!is.null(opt$infile)&&file.exists(opt$infile)){
  rdt <- read.table(opt$infile, stringsAsFactors = F, comment.char = '', header = T, check.names = F, sep = '\t')
}else{
  cat("WARNING: opt$infile not exist!\n\n")
  quit()
}

names(rdt) <- c('Class', names(rdt)[-1])
rdt$Class <- factor(rdt$Class, levels = unique(rdt$Class))

#step1 split dataset into trainning and testing files
set.seed(123)
yn <- sample(2, nrow(rdt), replace = TRUE, prob = c(opt$train_set_proportion, 1-opt$train_set_proportion))
train.dt <- rdt[yn == 1, ]
test.dt <- rdt[yn == 2, ]

###############################
###cross validation############
###finding best mtry and ntree
##############################
ntrees <- seq(100, 1000, by = 100)
err.cv <- vector()
#step2 determine the best choice of mtry and ntree using cross validation (time-consuming step)
for(i in 1:length(ntrees)){
  cv.tmp <- rfcv(train.dt[, -1], train.dt[,1], ntree = ntrees[i], cv.fold = opt$k_fold)
  err.cv <- cbind(err.cv, cv.tmp$error.cv)
}
err.cv.mat <- data.frame(rownames(err.cv), err.cv)

colnames(err.cv.mat) <- c('mtry.ntree', ntrees)

write.csv(err.cv.mat, file = paste0(opt$outdir, '/err.cv.csv'), quote = F, row.names = F)

colnames(err.cv.mat) <- c('mtry', ntrees)
m.ecm <- melt(err.cv.mat, id.vars = 'mtry', variable.name = 'ntree', value.name = 'error.rate')

m.ecm$mtry <- factor(m.ecm$mtry, levels = as.character(sort(as.numeric(levels(m.ecm$mtry)))))

min.v <- as.numeric(m.ecm[which.min(m.ecm$error.rate), ])
min.err <- min.v[3]
bmtry <- as.numeric(levels(m.ecm$mtry)[min.v[1]])
bntree <- as.numeric(levels(m.ecm$ntree)[min.v[2]])
cat("The best combination for mtry and ntree is", bmtry, 'and', bntree, ', error rate is ', min.err, ".\n")
p1 <- 
  ggplot(m.ecm, aes(mtry, ntree, size = error.rate))+
  geom_point(color = 'steelblue') + 
  annotate('text', label = paste0('Min error rate\n(', bmtry, ' ', bntree, ')', min.err), 
           x = min.v[1],  y = min.v[2],
           color = 'red', size = 3.5)+
  theme(legend.position = 'right')

ggsave(plot = p1, filename = paste0(opt$outdir, '/error.rate.cv.pdf'), height = 5, width = 7)

#step3 perform random forest using best parameters (mtry and ntree)
if(!is.null(opt$mtry) & !is.null(opt$ntree)){
  rf <- randomForest(train.dt[, -1], train.dt[,1], 
                     mtry = opt$mtry, ntree = opt$ntree,
                     importance = T, proximity = T )
}else{
  rf <- randomForest(train.dt[, -1], train.dt[,1], 
                     mtry = bmtry, ntree = bntree,
                     importance = T, proximity = T )
}

###out-of-bagging-error rate
oob.ed<- melt(rf$err.rate, varnames = c('Trees', 'Type'), value.name = "Error")
p2 <- 
  ggplot(oob.ed, aes(Trees, Error)) +
  geom_line(aes(color = Type), size = 1)
ggsave(plot = p2, filename = paste0(opt$outdir, '/oob.err.pdf'), height = 5, width = 7)

# hist(treesize(rf),
#      xlab = "No. of Variables for the Trees",
#      main = '', col = "green")

#step4 feature selection, MeanDecreaseAccuracy and MeanDecreaseGini
imp <- as.data.frame(importance(rf))
imp.out <- data.frame(Vars = rownames(imp), 
                      MeanDecreaseAccuracy = imp$MeanDecreaseAccuracy, 
                      MeanDecreaseGini = imp$MeanDecreaseGini)
write.table(imp.out, file = paste0(opt$outdir, '/important_vars.xls'), row.names = F, quote = F, sep = '\t')

imp.pd <- data.frame(Taxa = rownames(imp), 
                     MeanDecreaseAccuracy = imp$MeanDecreaseAccuracy, 
                     MeanDecreaseGini = imp$MeanDecreaseGini)
opt$top <- ifelse(opt$top >= nrow(imp.pd), nrow(imp.pd), opt$top)
temp <- imp.pd[order(-imp.pd$MeanDecreaseAccuracy), ][1:opt$top, ]
p3 <- 
  ggplot(temp, aes(reorder(Taxa, MeanDecreaseAccuracy), MeanDecreaseAccuracy)) +
  geom_bar(stat = 'identity', position = 'dodge',
           fill = 'steelblue', color = 'grey') + 
  xlab('') + 
  coord_flip()
ggsave(plot = p3, filename = paste0(opt$outdir, '/MeanDecreaseAccuracy.pdf'), height = 7, width = 6)

temp <- imp.pd[order(-imp.pd$MeanDecreaseGini), ][1:opt$top, ]
p4 <- 
  ggplot(temp, aes(reorder(Taxa, MeanDecreaseGini), MeanDecreaseGini)) +
  geom_bar(stat = 'identity', position = 'dodge',
           fill = 'steelblue', color = 'grey') + 
  xlab('') + 
  coord_flip()
ggsave(plot = p4, filename = paste0(opt$outdir, '/MeanDecreaseGini.pdf'), height = 7, width = 6)

#step5 MDS plot based on proximity similarity matrix
dm <- dist(1-rf$proximity)
mds.res <- cmdscale(dm, eig = TRUE, x.ret = TRUE)
egi.perc <- round(mds.res$eig/sum(mds.res$eig)*100, 2)[1:2]
mds1 <- paste0('MDS1 (', egi.perc[1], '%)')
mds2 <- paste0('MDS2 (', egi.perc[2], '%)')
pd <- data.frame(Sample = rownames(mds.res$points),
                       MDS1 = mds.res$points[,1],
                       MDS2 = mds.res$points[,2],
                       Level = train.dt$Class)
p5 <- ggplot(pd, aes(MDS1, MDS2, color = Level, fill = Level)) + 
  geom_point(aes(color = Level), size=5) +
  theme_bw() +
  xlab(mds1) +
  ylab(mds2) + 
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

ggsave(p5, filename = paste0(opt$outdir, '/MDS.rf.pdf'), height = 8, width = 8)
# MDSplot(rf, train.dt$Class, cex = 1.2, pch = 19)

#step6 prediction for trainning set and testing set based the rf model 
p1 <- predict(rf, train.dt)
p2 <- predict(rf, test.dt)
# cfm1 <- confusionMatrix(p1, train.dt[,1])
# cfm2 <- confusionMatrix(p2, test.dt[,1])
# acc1 <- sum(diag(cfm1$table))/sum((cfm1$table))
# acc2 <- sum(diag(cfm2$table))/sum((cfm2$table))
acc1 <- sum(p1 == train.dt[,1]) / length(train.dt[,1])
acc2 <- sum(p2 == test.dt[,1]) / length(test.dt[,1])
cat('The accuracy of this randomerorest model for trainning data:', acc1, "\n")
cat('The accuracy of this randomerorest model for testing data:', acc2, "\n")

pp1 <- predict(rf, train.dt, type = 'prob')
ppd1 <- data.frame(Class = as.vector(train.dt$Class), pp1)
m.ppd1 <- melt(ppd1, id.vars = 'Class', variable.name = 'Pred_Class', value.name = 'Prob')
m.ppd1$Class <- factor(m.ppd1$Class, levels = unique(rdt$Class))
m.ppd1$Pred_Class <- factor(m.ppd1$Pred_Class, levels = unique(rdt$Class))
m.ppd1$Consistence <- ifelse(m.ppd1$Class == m.ppd1$Pred_Class, 'Y', 'N')

bdt1 <- boxdata_get(ppd1)

pp2 <- predict(rf, test.dt, type = 'prob')
ppd2 <- data.frame(Class = as.vector(test.dt$Class), pp2)
m.ppd2 <- melt(ppd2, id.vars = 'Class', variable.name = 'Pred_Class', value.name = 'Prob')
m.ppd2$Class <- factor(m.ppd2$Class, levels = unique(rdt$Class))
m.ppd2$Pred_Class <- factor(m.ppd2$Pred_Class, levels = unique(rdt$Class))
m.ppd2$Consistence <- ifelse(m.ppd2$Class == m.ppd2$Pred_Class, 'Y', 'N')
bdt2 <- boxdata_get(ppd2)
write.table(ppd1, file = paste0(opt$outdir, '/trainset_prediction_prob.xls'), sep = '\t', quote = F, row.names = F)
write.table(ppd2, file = paste0(opt$outdir, '/testset_prediction_prob.xls'), sep = '\t', quote = F, row.names = F)

p6.1 <- 
  ggplot(bdt1, aes(Class, Prob)) +
  geom_boxplot(aes(fill = Class))
ggsave(p6.1, filename = paste0(opt$outdir, '/boxplot_prob_trainset1.pdf'), height = 6, width = 7)

p6.2 <- 
  ggplot(bdt2, aes(Class, Prob)) +
  geom_boxplot(aes(fill = Class))
ggsave(p6.2, filename = paste0(opt$outdir, '/boxplot_prob_testset1.pdf'), height = 6, width = 7)

p6.3 <- 
  ggplot(m.ppd1, aes(Class, Prob)) + 
  geom_boxplot()+
  facet_wrap(Pred_Class~.) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p6.3, filename = paste0(opt$outdir, '/boxplot_prob_trainset2.pdf'), height = 8, width = 8)

p6.4 <- 
  ggplot(m.ppd2, aes(Class, Prob)) + 
  geom_boxplot()+
  facet_wrap(Pred_Class~.) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p6.3, filename = paste0(opt$outdir, '/boxplot_prob_testset2.pdf'), height = 8, width = 8)

#step7 multi-group roc
###correct the error of auc 0.500
multiple_roc_plot(ppd2$Class, ppd2[,2], outdir = paste0(opt$outdir, '/ROC'))
cat('####################FINISHED:', date(),'####################\n')
###END###
