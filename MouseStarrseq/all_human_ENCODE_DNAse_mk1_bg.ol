; oligo-analysis  -v 2 -quick -1str -return freq,occ -i wgEncodeRegDnaseClusteredV2.fa -l 2 -noov -o all_DNAse.ol
; Citation: van Helden et al. (1998). J Mol Biol 281(5), 827-42. 
; Program version              	1.169
; Quick counting mode          
; Oligomer length              	2
; Input file                   	wgEncodeRegDnaseClusteredV2.fa
; Input format                 	fasta
; Output file                  	all_DNAse.ol
; Discard overlapping matches
; Counted on a single strand
; Background model             	Bernoulli
; Background estimation method 	
; Sequence type                	DNA
; Nb of sequences              	1281988
; Sum of sequence lengths      	387138495
; discarded residues           	NA (quick mode)	 (other letters than ACGT)
; discarded occurrences        	NA (quick mode)	 (contain discarded residues)
; nb possible positions        	NA (quick mode)
; total oligo occurrences      	385856225
; total overlapping occurrences	26888781
; total non overlapping occ    	358967444
; alphabet size                	4
; nb possible oligomers        	16
; oligomers tested for significance	0
;
; column headers
;	1	seq            	oligomer sequence
;	2	id             	oligomer identifier
;	3	obs_freq       	observed relative frequency
;	4	occ            	observed occurrences
;	5	ovl_occ        	number of overlapping occurrences (discarded from the count)
;	6	forbocc        	forbidden positions (to avoid self-overlap)
#seq	id	obs_freq	occ	ovl_occ	forbocc
aa	aa	0.0545041899998	21030781	7149412	21030781
ac	ac	0.0510507741582	19698259	0	19698259
ag	ag	0.0775606198915	29927248	0	29927248
at	at	0.0549440921939	21200520	0	21200520
ca	ca	0.0772630245890	29812419	0	29812419
cc	cc	0.0550307876982	21233972	6325105	21233972
cg	cg	0.0172122271709	6641445	0	6641445
ct	ct	0.0775392362790	29918997	0	29918997
ga	ga	0.0624745007029	24106175	0	24106175
gc	gc	0.0585621444879	22596568	0	22596568
gg	gg	0.0550163030284	21228383	6323177	21228383
gt	gt	0.0511378558166	19731860	0	19731860
ta	ta	0.0437937083949	16898075	0	16898075
tc	tc	0.0624039796170	24078964	0	24078964
tg	tg	0.0774273008036	29875806	0	29875806
tt	tt	0.0543932445304	20987972	7091087	20987972
; Host name	cn-r7-27
; Job started	2014-06-12.183458
; Job done	2014-06-12.183548
; Seconds	3.62
;	user	3.62
;	system	0.41
;	cuser	45.79
;	csystem	0.8
