# mkdir -p ./traindir_train_334_1700

rm(list = ls())

file.prefix.train="traindir_train_334_1700/"

unlink(file.prefix.train)

dir.create(file.prefix.train, showWarnings=FALSE)

feature_matrix=ledpred2::mapFeaturesToCRMs(URL='http://ifbprod.aitorgonzalezlab.org/map_features_to_crms.php', positive.bed='334level3.bed', negative.bed='1700level1.bed', pssm='../data/drosoVertebrateJaspar_flyFactor_hPDI.tf', background.freqs='2nt_upstream_Mus_musculus_EnsEMBL-noov-2str.freq', genome='mm9', my.values=list.files('peaks', full.names=TRUE), crm.feature.file='traindir_train_334_1700/crm_features_noscale.tab', stderr.log.file = 'traindir_train_334_1700/stderr.log', stdout.log.file = 'traindir_train_334_1700/stdout.log' )
dat=read.table('traindir_train_334_1700/crm_features_noscale.tab');
dat_mod=as.data.frame(apply(dat, 2, function(x) x/sqrt(sum(x^2))));
dat_mod[,1]=dat[,1];
write.table(dat_mod, file='traindir_train_334_1700/crm_features.tab', quote=F, sep='\t', col.names=NA)

crm.feature.file = file.path(file.prefix.train, 'crm_features.tab')
crms=read.table(crm.feature.file)
y=crms[,1]
x=crms[,-1]
valid.times = 5
step.nb=20
feature.nb.vector = seq(from = 10, to = (ncol(x) - 1), by = step.nb)
ledpred.obj=ledpred2::ledpred(x=x, y=y, valid.times=valid.times, feature.nb.vector=feature.nb.vector, file.prefix = file.prefix.train)

