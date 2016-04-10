rm(list = ls())

file.prefix.train="traindir_train_114_1140_test_29_290/"
unlink(file.prefix.train)
dir.create(file.prefix.train, showWarnings=FALSE)

feature_matrix=ledpred2::mapFeaturesToCRMs(URL='http://ifbprod.aitorgonzalezlab.org/map_features_to_crms.php', positive.bed='114CADenhancers.bed', negative.bed='1140CADnegative.bed', pssm='insect_570.tf', background.freqs='../data/2nt_intergenic_Drosophila_melanogaster.freq', genome='dm3', ngs=list.files('ngs', full.names=TRUE), crm.feature.file='traindir_train_114_1140_test_29_290/crm_features_noscale.tab', stderr.log.file = 'traindir_train_114_1140_test_29_290/stderr.log', stdout.log.file = 'traindir_train_114_1140_test_29_290/stdout.log' )

crm.feature.file = file.path(file.prefix.train, 'crm_features_noscale.tab')
crms=read.table(crm.feature.file)
y=crms[,1]
x=crms[,-1]
valid.times = 5
step.nb=20
feature.nb.vector = seq(from = 10, to = (ncol(x) - 1), by = step.nb)
ledpred.obj=ledpred2::ledpred(x=x, y=y, valid.times=valid.times, feature.nb.vector=feature.nb.vector, file.prefix = file.prefix.train)

file.prefix.test="testdir_train_114_1140_test_29_290/"
unlink(file.prefix.test)
dir.create(file.prefix.test, showWarnings=FALSE)

feature_matrix=ledpred2::mapFeaturesToCRMs(URL='http://ifbprod.aitorgonzalezlab.org/map_features_to_crms.php', positive.bed='29_290CAD.bed', genome='dm3', pssm='insect_570.tf', background.freqs='../data/2nt_intergenic_Drosophila_melanogaster.freq', ngs=list.files('ngs', full.names=TRUE), feature.ranking = read.table(file.path(file.prefix.train, '_feature_ranking.txt'), header=T), feature.nb = read.table(file.path(file.prefix.train, '_best.feature.nb.txt'))[1,1], crm.feature.file='testdir_train_114_1140_test_29_290/crm_features_noscale.tab', stdout.log.file = 'testdir_train_114_1140_test_29_290/stdout.log', stderr.log.file = 'testdir_train_114_1140_test_29_290/stderr.log')

crm.feature.file = file.path(file.prefix.test, 'crm_features_noscale.tab')
crms_test=read.table(crm.feature.file)
x=crms_test[,-1]

ledpred.obj <- get(load(file.path(file.prefix.train, '_ledpred.rda')))
model=ledpred.obj$model.obj
selected.features = as.character(ledpred.obj$feature.ranking$FeatureName[1:ledpred.obj$best.feature.nb])
x = x[,selected.features]

score.file=file.path(file.prefix.test, "scores.txt")
scores = ledpred2::scoreData(x, model=model, score.file = score.file)

# Convert scores.txt to scores.bed

score.file.bed=file.path(file.prefix.test, "scores.bed")
scores = read.table(score.file,header=F)
scores=cbind(data.frame(do.call('rbind', strsplit(as.character(scores[,1]),'_',fixed=TRUE))), scores$V2)
scores=scores[,-1]
scores$X5="."
scores=cbind(scores, data.frame(X6="+"))
write.table(scores, file=score.file.bed, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# Compute ROC curve using known positive in test data

positive.bed.fn = '29CADenhancers.bed'
negative.bed.fn = '290CADnegative.bed'
score.file.bed = 'scores.bed'

devtools::install_github('davetang/bedr')
require(bedr)
scores <- bed_to_granges(score.file.bed)

labels.pos <- bed_to_granges(positive.bed.fn)
labels.pos$label=1
labels.neg <- bed_to_granges(negative.bed.fn)
labels.neg$label=-1
labels=c(labels.pos, labels.neg)
scores.df=as.data.frame(scores)
labels.df=as.data.frame(labels)

roc.file=file.path(file.prefix.test, "tpr_fpr.png")
merged=merge(labels.df, scores.df, by.x=c('seqnames', 'start', 'end'), by.y=c('X2', 'X3', 'X4'))
require(ROCR)
merged$label <- as.factor(merged$label)
pred <- prediction(merged$score, merged$label, label.ordering = c(-1, 1))
perf <- performance(pred, measure = 'tpr', x.measure = 'fpr') 
auc = round(ROCR::performance(pred, 'auc')@y.values[[1]], digits=2);
png(roc.file)
plot(perf, col=rainbow(10))
title(paste('AUC: ', auc))
garb=dev.off()

