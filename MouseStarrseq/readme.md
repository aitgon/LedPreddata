# LedPred2 commands

To run the training on all 334 enhancers, use the commands in this file __traindir_train_334_1700.R__

To run the training on 267 enhancers, and test on remaining 67 enhancers, use the commands in this file __traindir_train_267_1360_test_67_340.R__


# Commands for ROC curve on test enhancers

Convert txt to bed with this awk command

~~~
awk 'BEGIN {FS="[_ ]"; OFS="\t"} {print $2,$3,$4,".",$NF,"+"}' "scores.txt" >scores.bed
~~~

And then run the commands here in the test folder: __test_prediction_roc.R__

# Data description

[Data ref](http://www.ncbi.nlm.nih.gov/pubmed/25872643)

3392 crms strong activity
5374 total

__Genome mm9__

# Original email from Salva

Salut Denis et Aitor,
Comme nous l’avons discuté je vous envoie un tableau de données avec le CapStarr-seq et des données épigénomiques.
 
Petit résumé de ce qu’il contient :
 
Uniquement des CRMs distaux qui overlappent avec au moins un TF et la DNaseI dans les thymocytes DP.
 
Le fold change du CApStarrseq (log2 ou pas) dans la lignée P5424 (proche des thymocytes DP)
 
Le classement en fonction de l’activité Capstarrseq (1=inactive ; 2=weak ; 3=strong)
 
Le signal ChIPseq pour chaque TF dans les Thymocytes DP
 
Le signal ChIPseq pour les modification d’histones et la Pol II dans la lignée P5424
 
 
Aurélien peux-tu nous dire si les signaux ChIP-seq sont en Log2 ou pas ?
 
 
Si vous avez des questions n’hésitez pas à nous demander.
 
A+
salva


