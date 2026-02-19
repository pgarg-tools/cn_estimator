library(tidyverse)
library(data.table)
library(reshape2)


genReference <- function(file){
  
  ### This function calculates mean and median per Loci to be considered as Reference Sample.
    data = read_rds(file) %>% 
      filter(sampleID %in% meta$sampleID) %>%
      filter(sampleID %in% qc$sampleID) %>%
      melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") 
    
    #            sampleID             Loci     CN
    # 1 1000250_23193_0_0 chr1:10000-15000 3.0462
    # 2 1000391_23193_0_0 chr1:10000-15000 3.2696
    # 3 1000712_23193_0_0 chr1:10000-15000 2.7744
    
    
    stat = data %>%
      group_by(Loci) %>%
      summarise(
        Mean = round(mean(CN),4),
        Median = round(median(CN), 4),
        Quant_0_01 = round(quantile(CN, probs = 0.0001),4),
        Quant_99_99 = round(quantile(CN, probs = 0.9999),4)
      ) %>% ungroup %>% as.data.frame
    
    stat_numsamples = data %>% inner_join(s, by  = "Loci") %>%
      group_by(Loci) %>%
      summarize(
        total = n(), # total number samples per Loci
        range = sum(CN >= Quant_0_01  & CN <=  Quant_99_99) # numSamples within quantiles per Loci
      ) %>% ungroup %>%
      as.data.frame
    
    
    stat =  stat %>% inner_join(stat_numsamples, by = "Loci")
    ref = stat %>% mutate(diff = Quant_99_99 - Quant_0_01 ) 
    
    #                   Loci   Mean Median Quant_0_01 Quant_99_99 total range   diff
    # 1     chr1:10000-15000 3.3226 3.2786     1.7925      7.8280 41615 41605 6.0355
    # 2   chr1:100000-105000 2.4815 2.4453     0.6629      5.2340 41615 41605 4.5711
    # 3 chr1:1000000-1005000 1.9981 1.9979     1.6887      2.9910 41615 41605 1.3023
    
    ref 
    
}

filterReference <- function(ref){
  ### This function filters Loci and keeps only invariant Loci per cohort
  ### Filtering criteria:
  ###    Diff between quantiles <0.9
  ###    Mean CN > =1.8 and <=2.2
  ###    Fraction of samples between quantiles  >99.9%
  ### This make sures that majority of the samples have relative CN near 2.
   
  ref = ref %>%
    separate(Loci, sep = ":|-", into = c("chr","start","end"), convert = T, remove = F)  %>% 
    filter(grepl("^chr[0-9]+$", chr)) %>% group_by(cohort) %>%
    filter(diff<0.9 & Mean >= 1.8 & Mean <= 2.2 & range/diff > 0.999) %>%
    ungroup %>% as.data.frame %>%
    group_by(Loci) %>%
    mutate(n = n()) %>%
    ungroup %>% as.data.frame
  
  ref %>% filter(n==5)  ### Make filtered Loci is invariant in all subbatches. In UK biobank I had 5 subcohorts
}

### store the Reference sample stats for all chromosomes in one file per subcohort. 

getCorFac <- function(cohort, chr){
  
  #### This cunction calculates correction factor for genomic waves
  #### Run this for each batch across all chromosomes. 
  print(cohort)
  
  band_filter = fread("Window_1kb.2mb.NoSatellite.NoSimpleRep.NoCentromere.txt", header = F, sep = "\t") %>% 
    as.data.frame %>% 
    rename(chr = 1, start = 2, end = 3, Loci = V4, Band = V5)
  
  ref_sample = REF_SAMPLE %>% rename(COHORT = cohort) %>% filter(COHORT == cohort) %>%
    select(Loci, Mean, Quant_0_01, Quant_99_99, COHORT) %>% rename(REF = Mean)
  
  
  print(head(ref_sample))
  
  files = list.files(paste0("../Data_5kb/",cohort), pattern = ".rds", full = T) %>% 
    grep(paste0(chr, "\\."), ., value = T)
  
  print(length(files))
  stat = mclapply(files, function(f){
    print(f)
    d = read_rds(f)  %>%
      melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") 
    
    coord = band_filter
    
    #### Since we selected invariant regions for reference samples close to 2, 
    #### we calculate genomic wave correction factor by calculating the mean deviation of sample CN from ref CN of 2
    d = d %>%
      inner_join(ref_sample, by = "Loci") %>%
      filter(CN >= Quant_0_01 & CN <=  Quant_99_99) %>%
      mutate(diff = CN - 2) %>%
      inner_join(coord, by = "Loci") %>%
      group_by(sampleID, Band) %>%
      summarize(
        MeanDiff = round(mean(diff),4),
        NumLoci = n(),
      ) %>%
      ungroup %>% as.data.frame
    d
  }, mc.cores = 8)
  
  
  ### Calculate genomic wave correction factor for each 2 mb band
  stat = stat %>% rbindlist %>% as.data.frame %>%
    group_by(sampleID, Band) %>%
    summarize(corfac = round(sum(MeanDiff * NumLoci)/sum(NumLoci),4), 
              NumLoci = sum(NumLoci) ) %>%
    ungroup %>% as.data.frame
  
  stat
  
  #stat
}


update_corfac <- function(x){
  
  ### this function updates the correction factors for those 2mb bands with missing values
  ### it includes the data for the two adujusant bands to calculate the correction factors the band with missing values.
  print(f)

  band = fread("Window_1kb.2mb.txt", header = F, sep= "\t") %>% as.data.frame %>%
    rename(chr = 1, start = 2, end = 3, Loci = V4, Band = V5)
  
  x = x %>% mutate(cohort = sub(".chr[0-9]+.perSample.CorrectionFactor.txt.gz","",basename(f) ))
  chrSelect = basename(f) %>% str_split(., "\\.") %>% map_chr(3)
  
  samples = x %>% select(sampleID) %>% distinct

  y = band %>%
    select(-Loci) %>%
    distinct %>%
    separate(Band, sep = ":|-", into = c("chr","start", "end"), convert = T, remove = F) %>% 
    filter(chr == chrSelect) %>% distinct %>%
    crossing(samples) %>%
    #  # filter(sampleID %in% select_samples) %>%
    arrange(chr, start, end) %>%
    group_by(sampleID, chr) %>%
    mutate(id = 1:n()) %>%
    left_join(x, by = c("Band", "sampleID")) %>% 
    mutate(corfac_numLoci = corfac * NumLoci) %>%
    mutate(code = ifelse(!is.na(NumLoci) & NumLoci >=50, "O", "P")) %>%
    mutate(lower_index = ifelse(code == "P", NA, id)) %>%
    mutate(upper_index = ifelse(code == "P", NA, id)) %>%
    ungroup %>% as.data.frame
  
  z = y %>%
    group_by(sampleID, chr) %>%
    arrange(start, end) %>%
    tidyr::fill(lower_index) %>%
    tidyr::fill(upper_index, .direction = "up") %>%
    mutate(lower_index = ifelse(code == "P" & id == 1 & is.na(lower_index), 1, lower_index)) %>%
    mutate(upper_index = ifelse(code == "P" & id == max(id) & is.na(upper_index), max(id), upper_index)) %>%
    tidyr::fill(lower_index) %>%
    tidyr::fill(upper_index, .direction = "up") %>%
    
    mutate(corfac2 = ifelse(id == lower_index & id == upper_index,
                            corfac,
                            map2_dbl(lower_index, upper_index, ~ sum(corfac_numLoci[.x:.y], na.rm = T)/sum(NumLoci[.x:.y], na.rm = T)) )
    ) %>% select(-corfac_numLoci) %>%
    # mutate(corfac2 = round(corfac2, 2)) %>%
    ungroup %>% as.data.frame() %>% mutate(corfac2 = round(corfac2, 4))
  z
}


getData <- function(file){
  cat(file); cat("\n")
  read_rds(file) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") %>% 
    mutate(cohort = basename(dirname(file)))
}


genomicWaveCorr <- function(file){
  
  data = getData(file) %>% rbindlist %>% as.data.frame  
  
  band = fread("Window_1kb.2mb.txt", header = F, sep= "\t") %>% as.data.frame %>%
    rename(chr = 1, start = 2, end = 3, Loci = V4, Band = V5) %>%
    select(Band, Loci) %>% filter(Loci %in% unique(data$Loci))
  
  coord = data %>% select(Loci) %>% distinct %>% 
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), remove = F, convert = T)
  
  CorrFac = getCorFac("SC.200K", "chr1") %>%  ### please modify it accordingly
    update_corfac %>%
    select(Band, sampleID, corfac2, cohort) %>%
    rename(corfac = corfac2) ### use updated genomic wave correction factor
  
  CORFAC = coord %>%
    inner_join(band, by = "Loci") %>%
    inner_join(CorrFac, by = c("Band")) %>%
    select(-cohort) %>% select(Loci, sampleID, corfac)
  
  old = nrow(data)
  data = data %>% 
    inner_join(CORFAC, by = c("sampleID","Loci"))
  new = nrow(data)
  
  if(new != old) { stop("ERROR: Number of rows in befor and after join do not match!!!") }
  data = data %>% mutate(CN = CN - corfac)
  data
}


batchCorrection <- function(data){
  stat_perCohort = data %>% 
    group_by(Loci, cohort) %>%
    summarise(median_cohort = round(median(CN), 4)) %>%
    ungroup %>% as.data.frame 
  
  stat_all = data %>% group_by(Loci) %>%
    summarise(median_all = round(median(CN), 4)) %>%
    ungroup %>% as.data.frame
  
  stat = stat_perCohort %>% inner_join(stat_all, by = "Loci")
  
  data = data %>% 
    inner_join(stat, by = c("Loci", "cohort")) %>%
    mutate(CN2 = CN + median_all - median_cohort ) %>%
    mutate(CN2 = round(CN2, 4)) %>%
    select(Loci, sampleID, CN2 ) %>%
    spread(Loci, CN2)
  
  data
}


