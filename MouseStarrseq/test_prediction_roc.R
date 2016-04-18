positive.bed.fn = '../67_level3_starrseq.bed'
negative.bed.fn = '../340_level1_starrseq.bed'
scores.fn = 'scores.bed'
devtools::install_github('davetang/bedr')
require(bedr)
scores <- bed_to_granges(scores.fn)

labels.pos <- bed_to_granges(positive.bed.fn)
labels.pos$label=1
labels.neg <- bed_to_granges(negative.bed.fn)
labels.neg$label=-1
labels=c(labels.pos, labels.neg)
scores.df=as.data.frame(scores)
labels.df=as.data.frame(labels)

merged=merge(labels.df, scores.df, by=c('seqnames', 'start', 'end'))
#require(ROCR)
merged$label <- as.factor(merged$label)
pred <- ROCR::prediction(merged$score, merged$label, label.ordering = c(-1, 1))
perf <- ROCR::performance(pred, measure = 'tpr', x.measure = 'fpr') 
auc = round(ROCR::performance(pred, 'auc')@y.values[[1]], digits=2);
png('tpr_fpr.png')
plot(perf, col=rainbow(10))
title(paste('AUC: ', auc))
garb=dev.off()

