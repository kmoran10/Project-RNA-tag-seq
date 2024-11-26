

# Overlap analysis 



library(tidyverse)

LH <- readRDS("results/LH_limma_results.RDS")
DMH <- readRDS("results/DMH_limma_results.RDS")
ARC <- readRDS("results/ARC_limma_results.RDS")

LH_sigs <- LH %>% 
  filter(logFC >= .2 | logFC <= -.2) %>%
  filter(P.Value < 0.05) %>% 
  arrange(desc(abs(logFC))) %>% 
  select(symbol, logFC, P.Value, description) %>% 
  mutate(direction = if_else(logFC > 0, "Down", "Up"))

DMH_sigs <- DMH %>% 
  filter(logFC >= .2 | logFC <= -.2) %>%
  filter(P.Value < 0.05) %>% 
  arrange(desc(abs(logFC))) %>% 
  select(symbol, logFC, P.Value, description) %>% 
  mutate(direction = if_else(logFC > 0, "Down", "Up"))

ARC_sigs <- ARC %>% 
  filter(logFC >= .2 | logFC <= -.2) %>%
  filter(P.Value < 0.05) %>% 
  arrange(desc(abs(logFC))) %>% 
  select(symbol, logFC, P.Value, description) %>% 
  mutate(direction = if_else(logFC > 0, "Down", "Up"))


combined <- bind_rows(
  LH_sigs %>% mutate(source = "LH"),
  DMH_sigs %>% mutate(source = "DMH"),
  ARC_sigs %>% mutate(source = "ARC")
)



overlap_table <- combined %>%
  group_by(symbol, direction) %>%
  summarize(
    count = n_distinct(source),
    sources = paste(unique(source), collapse = ", "),
    description = paste(unique(description), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(count > 1) %>% # Remove cases with no overlap
  arrange(desc(count)) %>% 
  arrange(desc(direction)) %>% 
  select(symbol, direction, sources, description)

write.csv(overlap_table,"results/results_tables/overlap_table.csv", row.names = F)





conflict_table <- combined %>%
  group_by(symbol) %>%
  filter(n_distinct(direction) > 1) %>% # Keep symbols with conflicting directions
  summarize(
    sources_up = paste(source[direction == "Up"], collapse = ", "),
    sources_down = paste(source[direction == "Down"], collapse = ", "),
    description = paste(unique(description), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(sources_up != "" & sources_down != "") %>% # Ensure there's at least one source for both "Up" and "Down"
  arrange(desc(sources_up))

write.csv(conflict_table,"results/results_tables/conflict_table.csv", row.names = F)








