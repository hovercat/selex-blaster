setwd("/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all")

library(vroom)

#temp <- list.files(pattern="*.tsv", recursive = TRUE)
#motif_search_fimo <- mapply(vroom, temp, comment='#')
#cluster_names <- mapply(substr, temp, 0, 7)
#names(motif_search_fimo) <- cluster_names

selex_data <- read.delim("selex.aptamers.raw.csv")
selex_data_rounds <- selex_data %>%
  pivot_longer(-sequence, 
               names_to="round_id",
               values_to="seq_count") %>%
  group_by(round_id) %>%
  summarize(
    round_seq_count = sum(seq_count)
  )


library(dplyr)
library(tidyr)
tag_file <- read.delim("motif_analysis.tags.csv", sep="\t") %>%
  group_by(cluster_id) %>%
  arrange(tag)

#motif_search_fimo <- read.delim("./motif_analysis.csv", sep=",")
#motif_search_fimo_long <- motif_search_fimo %>%
#  pivot_longer(c(-cluster, -motif_id), names_to="round", values_to="perc")
motif_search_fimo_long <- vroom("./mast_cat.csv")
motif_search_fimo_long$X <- NULL


clusters <- motif_search_fimo_long %>% 
  inner_join(selex_data, by=c("sequence_name"="sequence")) %>%
  group_by(cluster) %>%
  summarize_at(
    vars(-sequence_name, -combined_p_val, -ev), sum, na.rm=TRUE
  )

clusters <- clusters %>%
  pivot_longer(-cluster, 
               names_to="round_id",
               values_to="seq_count")

clusters <- clusters %>% inner_join(selex_data_rounds, by=c("round_id"="round_id")) %>% mutate(seq_count_p = seq_count/round_seq_count, round_seq_count = NULL)


rounds = clusters %>%
  inner_join(selex_data_rounds) %>%
  select(round_id, round_seq_count) %>%
  unique()

rounds_names <- rounds$round_id %>% unique()
rounds_joiner <- data.frame(round_id = rounds_names[1:length(rounds_names)-1], next_round_id = rounds_names[2:length(rounds_names)])
rounds_joiner <- rounds_joiner %>% 
  inner_join(rounds) %>%
  inner_join(rounds, by=c("next_round_id"="round_id"))

clusters_joined <- 
  clusters %>%
  inner_join(rounds_joiner) %>%
  inner_join(clusters, by=c("cluster", "next_round_id" = "round_id")) %>%
  filter(
    seq_count.x > 10,
    seq_count.y > 10
  ) %>%
  mutate(
    enrichment = seq_count_p.y / seq_count_p.x,
    log_enrichment = log(enrichment, base=2)
  )

cluster_joined_centered <- clusters_joined %>%
  group_by(round_id) %>%
  mutate(
#    enrichment = enrichment / mean(enrichment), # centered
#    log_enrichment = log(enrichment, base=2) # centered
  )

mast_analysis <- cluster_joined_centered %>%
  group_by(cluster, round_id) %>%
  mutate(
    positive_round = next_round_id %in% c("R5", "R9"),
    negative_round = !positive_round,
    score = positive_round * log_enrichment  - negative_round * log_enrichment
  ) %>%
  group_by(cluster) %>%
  summarize(
    score_sum = sum(score),
    score_sd = sd(score),
    rounds=n()
  ) #%>% left_join(tag_file, by=c("cluster"="cluster_id")) %>%
  group_by(tag) %>%
  mutate(
    tag_in_clusters = n()
  ) %>%
  group_by(cluster) %>%
  mutate(
    tags_in_cluster = sum(!is.na(tag))
  ) %>%
  ungroup()


fisher <- function(a,b,c,d, alternative){
  data <- matrix(c(a,b,c,d),ncol=2)
  f <- fisher.test(data, alternative=alternative)
  c(f = f$p.value,
    OR = f$estimate)
}

mast_analysis_fisher <- cluster_joined_centered %>%
  group_by(cluster, round_id) %>%
  #filter(cluster == "c000000") %>%
  rowwise() %>%
  mutate(
    fi_p_gt = fisher(round_seq_count.x, seq_count.x, round_seq_count.y, seq_count.y, "greater")[[1]],
    fi_OT_gt = fisher(round_seq_count.x, seq_count.x, round_seq_count.y, seq_count.y, "greater")[[2]],
    fi_p_lt = fisher(round_seq_count.x, seq_count.x, round_seq_count.y, seq_count.y, "less")[[1]],
    fi_OT_lt = fisher(round_seq_count.x, seq_count.x, round_seq_count.y, seq_count.y, "less")[[2]]
  )

mast_analysis_fisher_mut <- mast_analysis_fisher %>%
 # filter(round_id == "R0") %>%
  mutate(
    positive_round = next_round_id %in% c("R4", "R8"),
    negative_round = !positive_round,
    score_p = positive_round * fi_p_gt,
    score_n = negative_round * fi_p_lt,
    score = score_n
    #score = log10(1 - score_p + score_n),
    #score =  (score_p > 0.9) + (score_p > 0.9999) + (score_p > (1- 10^(-20))) + (score_n > 0.99)
  ) 
mast_analysis_fisher_mut <- mast_analysis_fisher_mut %>%
  group_by(cluster) %>%
  summarize(
    score_mean = mean(score),
    score_sum = sum(score),
    score_prod = prod(score),
    score_sd = sd(score),
    rounds=n()
  ) 
#%>% left_join(tag_file, by=c("cluster"="cluster_id")) #%>%
  group_by(tag) %>%
  mutate(
    tag_in_clusters = n()
  ) %>%
  group_by(cluster) %>%
  mutate(
    tags_in_cluster = sum(!is.na(tag))
  ) %>%
  ungroup()

mast_analysis_summed <- mast_analysis_fisher_mut %>% inner_join(tag_file, by=c("cluster"="cluster_id")) %>% group_by(tag) %>%
  summarise(
    min_score = min(score_sum),
    mean_score = mean(score_sum),
    max_score = max(score_sum),
    mean_score_sd = sd(score_sum),
    score_sd_sd_min = min(score_sd),
    score_sd_sd_mean = mean(score_sd),
    score_sd_sd_max = max(score_sd),
    
    nr_of_clusters = n()
  )
#View(mast_analysis_summed)
tags <- paste0("EF", 701:718)
mas_tagged <- mast_analysis_fisher_mut # %>% filter(tag %in% tags)
library(ggplot2)

tagged_selex_data <- selex_data %>% inner_join(tag_file, by=c("sequence"="seq")) %>% #%>% mutate(cluster_id = NULL) %>%
  unique() %>% 
  filter(tag %in% tags) %>%
  group_by(tag)#%>%
  #slice_min(combined_p_val,n=5)
  #slice_min(ev,n=10)
 #slice_max(,n=10)

mas_selex <- tagged_selex_data %>% inner_join(mas_tagged, by=c("cluster_id"="cluster"))# %>%
  slice_sample(n=100)
  slice_min(score_sd, n=5) %>%
  #slice_max(score_mean, n=5)


mas_selex_long <- mas_selex %>%
  pivot_longer(
    c(-sequence, -tag, -cluster_id, -score_mean,-score_sum, -score_sd, -score_prod, -rounds, -ev, -combined_p_val),
    names_to="round",
    values_to="seq_count",
    names_prefix="R"
  ) %>%
  group_by(tag) %>%
  slice_max(score_mean, n=1)

gg1 <- ggplot(data=mas_selex, aes(x=tag)) +
  #geom_line(aes(y=score_sum), stat="identity") +
  geom_boxplot(aes(y=score_sum, group=tag, color="score_meanm")) +
  theme(axis.text.x = element_text(angle = 40), legend.position = c(0.9, 0.2)) #+
  #scale_y_log10()
  

gg1
  
