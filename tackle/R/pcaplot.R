library(tidyverse)
library(ggfortify)

## library(cluster)

exprs_long <- data %>% pivot_longer(c(-GeneSymbol, -GeneID)) %>%
  left_join(col_data, by='name')

forpca <- exprs_long %>% pivot_wider(id_cols=c(name, model, subtype), names_from=GeneID, values_from=value)
pca_res <- prcomp(forpca%>%select(-name, -model, -subtype), scale. = FALSE, center=TRUE)

p <- autoplot(pca_res, data=forpca, colour='model', shape='subtype', label=TRUE, label.repel=TRUE,
              label.label='name',
         frame=TRUE, frame.type='norm',
         label.size=3,
         size=2
         ) +
  ggplot2::theme_classic()
ggsave('abbvie_pca.pdf', p, width=7, height=7)
