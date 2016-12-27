#!/usr/bin/env Rscript

library(data.table)
library(dtplyr)
library(dplyr)
library(readr)
library(tidyr)

cp = data.table(read_tsv('~/data/pfitmap-eval/classified_proteins.prop_matching_ge_0.9.tsv'))

NrdGRE = cp %>% filter(psuperfamily=='NrdGRE')
org = cp %>% select(tdomain:tstrain) %>% distinct()

NrdGRE.domclass = NrdGRE %>% group_by(tdomain, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
write_tsv(NrdGRE.domclass, 'NrdGRE.domclass.tsv')
write_tsv(NrdGRE.domclass %>% spread(pclass, n_proteins), 'NrdGRE.domclass.wide.tsv')

NrdGRE.domsubclass = NrdGRE %>% group_by(tdomain, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
write_tsv(NrdGRE.domsubclass, 'NrdGRE.domsubclass.tsv')

NrdGRE.orderclass = NrdGRE %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
write_tsv(NrdGRE.orderclass, 'NrdGRE.orderclass.tsv')

NrdGRE.ordersubclass = NrdGRE %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
write_tsv(NrdGRE.ordersubclass, 'NrdGRE.ordersubclass.tsv')

NrdGREfer = cp %>% filter(psuperfamily %in% c('NrdGRE', 'Ferritin-like'))

NrdGREfer.domclass = NrdGREfer %>% group_by(tdomain, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
write_tsv(NrdGREfer.domclass, 'NrdGREfer.domclass.tsv')

NrdGREfer.domsubclass = NrdGREfer %>% group_by(tdomain, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
write_tsv(NrdGREfer.domsubclass, 'NrdGREfer.domsubclass.tsv')

NrdGREfer.orderclass = NrdGREfer %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
write_tsv(NrdGREfer.orderclass, 'NrdGREfer.orderclass.tsv')

NrdGREfer.ordersubclass = NrdGREfer %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
write_tsv(NrdGREfer.ordersubclass, 'NrdGREfer.ordersubclass.tsv')