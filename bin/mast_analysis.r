setwd("/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all")

library(vroom)

temp <- list.files(pattern="*.tsv", recursive = TRUE)
motif_search_fimo <- mapply(vroom, temp, comment='#')
cluster_names <- mapply(substr, temp, 0, 7)
names(motif_search_fimo) <- cluster_names

selex_data <- vroom("../selex.aptamers.raw.csv")

motif_search_fimo$c000000$
  vroom


library(dplyr)
library(tidyr)
tag_file <- read.delim("motif_analysis.tags.csv", sep="\t") %>%
  group_by(cluster_id, motif_id) %>%
  arrange(tag)

#motif_search_fimo <- read.delim("./motif_analysis.csv", sep=",")
#motif_search_fimo_long <- motif_search_fimo %>%
#  pivot_longer(c(-cluster, -motif_id), names_to="round", values_to="perc")
motif_search_fimo_long <- read.delim("./motif_analysis.long.csv", sep=",")

rounds = motif_search_fimo_long$round_id %>% unique()
rounds_joiner = data.frame(round_id = rounds[1:length(rounds)-1], next_round_id = rounds[2:length(rounds)])


msfl_enrichment_roundwise <- motif_search_fimo_long %>% 
  inner_join(rounds_joiner) %>%
  inner_join(motif_search_fimo_long, by=c("cluster", "motif_id", "next_round_id"="round_id"))  %>%
  mutate(
    enrichment = count_p.y / count_p.x,
    log_enrichment = log(enrichment, base=2)
  ) %>%
  group_by(round_id) %>%
  summarise(
    enrichment = mean(enrichment),
    log_enrichment = mean(log_enrichment)
  )

msfl_enrichment <- motif_search_fimo_long %>% 
  inner_join(rounds_joiner) %>%
  inner_join(motif_search_fimo_long, by=c("cluster", "motif_id", "next_round_id" = "round_id")) %>%
  mutate(
    enrichment = count_p.y / count_p.x,
    log_enrichment = log(enrichment, base=2)
  ) %>%
  group_by(round_id) %>%
  mutate(
    enrichment_centered = enrichment / mean(enrichment),
    log_enrichment_centered = log(enrichment_centered, base=2)
  )

msfl_rs <- msfl_enrichment %>% group_by(cluster, motif_id, round) %>%
  mutate(plus_round = next_round %in% c("X9", "X5"),
         neg_round = !plus_round,
         round_score = plus_round * log_enrichment_centered - neg_round * log_enrichment_centered)  %>%
  group_by(cluster, motif_id) %>%
  summarize(round_score_sum = sum(round_score), round_score_sd = sd(round_score)) %>%
  left_join(tag_file)
#View(msfl_rs)


motif_search_fimo_long %>% 
  inner_join(rounds_joiner) %>%
  inner_join(motif_search_fimo_long, by=c("cluster", "motif_id", "next_round_id" = "round_id"))


msfl_enrichment <- motif_search_fimo_long %>% 
  inner_join(rounds_joiner) %>%
  inner_join(motif_search_fimo_long, by=c("cluster", "motif_id", "next_round_id" = "round_id")) %>%
  filter(
    count_total.x > 25,
    count_total.y > 25
  ) %>%
  mutate(
    enrichment = count_p.y / count_p.x,
    #enrichment = (count_total.y / n.y) / (count_total.x / n.x),
    log_enrichment = log(enrichment, base=2)
  ) %>%
  group_by(round_id) %>%
  mutate(
    enrichment = enrichment / mean(enrichment), # centered
    log_enrichment = log(enrichment, base=2) # centered
  ) %>% 
  group_by(cluster, motif_id, round_id) %>%
  mutate(plus_round = next_round_id %in% c(5, 9),
         neg_round = !plus_round,
         round_score = plus_round * log_enrichment + (-1 * neg_round) * log_enrichment
         
         ) %>%
  ungroup() %>%
  group_by(motif_id, cluster) %>%
  summarize(
    round_score_sum = sum(round_score, na.rm=TRUE), 
    round_score_sd = sd(round_score, na.rm=TRUE),
    rounds=n(), 
    med_fimo_score=first(med_fimo_score.x), 
    sd_fimo_score=first(sd_fimo_score.x)
    
  ) %>%
  left_join(tag_file)# %>%
  group_by(tag) %>%
  summarize(
    mean_fimo_score = mean(med_fimo_score),
    sd_fimo_score = sd(med_fimo_score),
    sd_score = sd(round_score_sum),
    mean_score = mean(round_score_sum),
    max_score = max(round_score_sum),
    n=n()
  )
View(msfl_enrichment)
