##############################
###table with all variables###
##############################

library(dplyr)
library(reshape2)

updated_207_quant %>% 
        melt %>%
        group_by(variable) %>% 
        summarize(mean = mean(value,na.rm = TRUE)
                  , median = median(value,na.rm = TRUE)
                  , min = min(value,na.rm = TRUE)
                  , max = max(value,na.rm = TRUE)
                  ) %>% 
        write.csv(.,file = "~/R_Annie/final graphics/variable_overview.csv")

##############################
#####boxplots by cluster#####
##############################

#export as 25x25 inch pdf

library(ggplot2)
library(ggpubr)
library(tidyr)

#create comparisons
my_comparisons <- list( c("Huia/Whatipu", "Piha/Karekare"), c("Piha/Karekare", "North"), c("Huia/Whatipu", "North") )

#create boxplots
`updated_207` %>% 
  pivot_longer(-cluster) %>%
  ggplot(aes(factor(cluster), value, na.rm = TRUE)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(vars(name), scales = "free_y") +
  xlab("") + ylab("") +
  stat_compare_means(method='kruskal'
                     , size=3
                     , fontface='bold'
                     , color="Red"
                     , label.y.npc = "bottom"
                     , label.x.npc = "left"
                     , na.rm = TRUE
                     , aes(label = paste0("overall : ", ..p.signif.., "(", ..p.format.., ")"))) +
  stat_compare_means(test='kruskal'
                     , comparisons = my_comparisons
                     , size=3
                     , na.rm = TRUE
                     , aes(label = ..p.signif..)
                     , hide.ns=TRUE
                     )

##############################
######correlation matrix######
##############################

library(corrplot)

#overall significance test
sig <- corrplot::cor.mtest(updated_207_quant, conf.level = .95)

#overall Spearman correlation coefficients matrix
corrplot(cor(na.omit(updated_207_quant, method = "spearman"))
         , addCoef.col = 1
         , number.cex = 0.18
         , tl.cex = 0.25
         , cl.cex = 0.3
         , title="All Quantitative Variables Correlation Matrix"
         , cex.main=.8
         , type = 'lower'
         , mar=c(0,0,1,0) ##fixes margins
         , p.mat = sig$p
         , sig.level = .05 
         , insig = "blank"
)

#clusters significance tests 
sig_north <- corrplot::cor.mtest(updated_207_n, conf.level = .95)
sig_huia <- corrplot::cor.mtest(updated_207_hw, conf.level = .95)
sig_piha <- corrplot::cor.mtest(updated_207_pk, conf.level = .95)

#North cluster Spearman correlation coefficients matrix
library(corrplot)
corrplot(cor(na.omit(updated_207_n, method = "spearman"))
         , addCoef.col = 1
         , number.cex = 0.18
         , tl.cex = 0.25
         , cl.cex = 0.3
         , title="North Cluster Quantitative Variables Correlation Matrix"
         , cex.main=.5
         , type = 'lower'
         , mar=c(0,0,1,0) ##fixes margins
         , p.mat = sig_north$p
         , sig.level = .05 
         , insig = "blank"
)

#Huia/Whatipu cluster Spearman correlation coefficients matrix
library(corrplot)
corrplot(cor(na.omit(updated_207_hw, method = "spearman"))
         , addCoef.col = 1
         , number.cex = 0.18
         , tl.cex = 0.25
         , cl.cex = 0.3
         , title="Huia/Whatipu Cluster Quantitative Variables Correlation Matrix"
         , cex.main=.5
         , type = 'lower'
         , mar=c(0,0,1,0) ##fixes margins
         , p.mat = sig_huia$p
         , sig.level = .05 
         , insig = "blank"
)

#Piha/Karekare cluster Spearman correlation coefficients matrix
library(corrplot)
corrplot(cor(na.omit(updated_207_pk, method = "spearman"))
         , addCoef.col = 1
         , number.cex = 0.18
         , tl.cex = 0.25
         , cl.cex = 0.3
         , title="Piha/Karekare Cluster Quantitative Variables Correlation Matrix"
         , cex.main=.5
         , type = 'lower'
         , mar=c(0,0,1,0) ##fixes margins
         , p.mat = sig_piha$p
         , sig.level = .05 
         , insig = "blank"
)
