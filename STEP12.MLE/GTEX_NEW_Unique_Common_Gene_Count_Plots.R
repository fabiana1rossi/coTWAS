# Set working directory to project root
# setwd("path/to/project/root")

load("./predictions/LIBD_training/GTeX.v9/Combined/MLE_common.RData")
regions = names(Common.Predictions.logTPM)
combos = names(Common.Predictions.logTPM$dlpfc)

common.genes.plot = sapply(regions, function(region){
  combo.df = sapply(combos, function(combo){
    genes = names(Common.Predictions.logTPM[[region]][[combo]])
    gene.df = sapply(genes, function(gene){
      if(combo %in% c("EMI","EM",
                      "EI","MIC",
                      "MC","IC")){

      trans.on.cov.p = Common.Predictions.logTPM[[region]][[combo]][[gene]]$trans.on.cov.anova.p
      trans.on.cis.p = Common.Predictions.logTPM[[region]][[combo]][[gene]]$trans.on.cis.anova.p
      cis.on.cov.p   = Common.Predictions.logTPM[[region]][[combo]][[gene]]$cis.on.cov.anova.p
      trans.R        = Common.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.trans
      cis.R          = Common.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.cis
      cov.R          = Common.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.cov
      
      df = data.frame(
        trans.on.cov.p = trans.on.cov.p,
        trans.on.cis.p = trans.on.cis.p,
        cis.on.cov.p = cis.on.cov.p,
        trans.R = trans.R,
        cis.R = cis.R,
        cov.R = cov.R,
        trans.cis.R = trans.R - cis.R,
        trans.cov.R = trans.R - cov.R,
        cis.cov.R = cis.R - cov.R,
        combo = combo,
        region = region,
        gene = gene
      )
      }else{
        trans.on.cov.p = Common.Predictions.logTPM[[region]][[combo]][[gene]]$trans.on.cov.anova.p
        trans.on.cis.p = NA
        cis.on.cov.p = NA
        trans.R = Common.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.trans
        cis.R = NA
        cov.R = Common.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.cov
        
        df = data.frame(
          trans.on.cov.p = trans.on.cov.p,
          trans.on.cis.p = NA,
          cis.on.cov.p = NA,
          trans.R = trans.R,
          cis.R = NA,
          cov.R = cov.R,
          trans.cis.R = NA,
          trans.cov.R = trans.R - cov.R,
          cis.cov.R = NA,
          combo = combo,
          region = region,
          gene = gene
        ) 
      }
      
      
      },simplify = F, USE.NAMES =T)
      do.call(rbind,gene.df)
  },simplify = F, USE.NAMES =T)
    do.call(rbind,combo.df)
},simplify = F, USE.NAMES =T)
  
common.genes.plot.df = do.call(rbind,common.genes.plot)
  
common.genes.plot.df$Trans.on.cov.significant = ifelse(common.genes.plot.df$trans.on.cov.p < 0.05,
                                                       "Yes", "No")
common.genes.plot.df$Trans.on.cis.significant = ifelse(common.genes.plot.df$trans.on.cis.p < 0.05,
                                                       "Yes", "No")
common.genes.plot.df$Cis.on.cov.significant = ifelse(common.genes.plot.df$cis.on.cov.p < 0.05,
                                                       "Yes", "No")

common.genes.plot.df$Trans.on.cov.winner = ifelse(
  common.genes.plot.df$Trans.on.cov.significant == "Yes" & 
    common.genes.plot.df$trans.cov.R > 0 & ## Trans R higher than Cov R
    common.genes.plot.df$trans.R > 0, "Yes", "No") ## Trans R not negative

common.genes.plot.df$Trans.on.cis.winner = ifelse(
  common.genes.plot.df$Trans.on.cis.significant == "Yes" & 
    common.genes.plot.df$trans.cis.R > 0 & ## Trans R higher than Cis R
    common.genes.plot.df$trans.R > 0, "Yes", "No") ## Trans R not negative

common.genes.plot.df$Cis.on.cov.winner = ifelse(
  common.genes.plot.df$Cis.on.cov.significant == "Yes" & 
    common.genes.plot.df$cis.cov.R > 0 & # Cis R higher than Cis R
    common.genes.plot.df$cis.R > 0, "Yes", "No")## Cis R not negative

common.genes.plot.df$Best_Model <- ifelse(common.genes.plot.df$Trans.on.cov.significant == "Yes", "Trans",
               ifelse(common.genes.plot.df$Cis.on.cov.significant == "Yes" & 
                        common.genes.plot.df$Trans.on.cov.significant != "Yes", "Cis", "Cov"))

common.genes.plot.df$Best_Model_Details = ifelse(common.genes.plot.df$Trans.on.cov.significant == "Yes" & common.genes.plot.df$Trans.on.cis.significant == "Yes", "Trans_on_Cis_and_Cov",
               ifelse(common.genes.plot.df$Trans.on.cov.significant == "Yes" & common.genes.plot.df$Trans.on.cis.significant != "Yes", "Trans_on_Cov_only",
                      ifelse(common.genes.plot.df$Trans.on.cov.significant != "Yes" & common.genes.plot.df$Trans.on.cis.significant == "Yes", "Trans_on_Cis_only",
                             ifelse(common.genes.plot.df$Trans.on.cis.significant != "Yes" & common.genes.plot.df$Cis.on.cov.significant == "Yes", "Cis_on_Trans_and_Cov", 
                                    ifelse(common.genes.plot.df$Cis.on.cov.significant == "Yes", "Cis_on_Cov_only", "not_significant")))))

plot_data <- common.genes.plot.df %>%
  group_by(region, combo, Best_Model) %>%
  summarise(
    Number.of.Genes = n(),
    .groups = "drop"
  )
c = ggplot(plot_data, aes(y = Number.of.Genes, x = region, 
                         group= Best_Model, fill = Best_Model)) +
  geom_bar(stat = 'identity', position = "stack") + theme_bw() + 
  scale_fill_manual(values = c("cadetblue2", "plum2", "lightsalmon"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
svg("TWAS/Common_Genes_count.svg", width = 10, height = 6)
c
dev.off()

plot_data2 <- common.genes.plot.df %>%
  group_by(region, combo, Best_Model_Details) %>%
  summarise(
    Number.of.Genes = n(),
    .groups = "drop"
  )
c2 = ggplot(plot_data2, aes(y = Number.of.Genes, x = region, 
                          group= Best_Model_Details, fill = Best_Model_Details)) +
  geom_bar(stat = 'identity', position = "stack") + theme_bw() + 
  #scale_fill_manual(values = c("cadetblue2", "plum2", "lightsalmon"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

svg("TWAS/Common_Genes_count_details.svg", width = 10, height = 6)
c2
dev.off()


common.genes.plot.df$Model_Combo = paste0(common.genes.plot.df$combo, "_", 
                                           common.genes.plot.df$Best_Model)
common.genes.plot.df$Model_Combo = gsub("common_prediction_","",common.genes.plot.df$Model_Combo)
plot_data3 <- common.genes.plot.df %>%
  group_by(region, Model_Combo) %>%
  summarise(
    Number.of.Genes = n(),
    .groups = "drop"
  )
c3 = ggplot(plot_data3, aes(y = Number.of.Genes, x = region, 
                            group= Model_Combo, fill = Model_Combo)) +
  geom_bar(stat = 'identity', position = "stack") + theme_bw() + 
  # scale_fill_manual(values = c("cadetblue2", "plum2", "lightsalmon"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

svg("TWAS/Common_Genes_count_All_Combos.svg", width = 10, height = 6)
c3
dev.off()




load("./predictions/LIBD_training/GTeX.v9/Combined/MLE_unique.RData")
regions = names(Unique.Predictions.logTPM)
combos = names(Unique.Predictions.logTPM$dlpfc)
unique.genes.plot = sapply(regions, function(region){
  combo.df = sapply(combos, function(combo){
    genes = names(Unique.Predictions.logTPM[[region]][[combo]])
    gene.df = sapply(genes, function(gene){
        print(paste0(region, "_", combo))
        model.on.cov.p = Unique.Predictions.logTPM[[region]][[combo]][[gene]]$anova.p
        model.R = Unique.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.model
        cov.R = Unique.Predictions.logTPM[[region]][[combo]][[gene]]$adj.r.cov
        
        df = data.frame(
          Anova.p = model.on.cov.p,
          model.R = model.R,
          cov.R = cov.R,
          Rdiff = model.R - cov.R,
          combo = combo,
          region = region,
          gene = gene
        )
      
    },simplify = F, USE.NAMES =T)
    do.call(rbind,gene.df)
  },simplify = F, USE.NAMES =T)
  do.call(rbind,combo.df)
},simplify = F, USE.NAMES =T)

unique.genes.plot.df = do.call(rbind,unique.genes.plot)


unique.genes.plot.df$Model.on.cov.significant = ifelse(unique.genes.plot.df$Anova.p < 0.05,
                                                       "Yes", "No")
unique.genes.plot.df$Model.on.cov.winner = ifelse(
  unique.genes.plot.df$Model.on.cov.significant == "Yes" & 
    unique.genes.plot.df$Rdiff > 0 & 
    unique.genes.plot.df$model.R > 0, "Yes", "No")

plot_data <- unique.genes.plot.df %>%
  group_by(region, combo, Model.on.cov.winner) %>%
  summarise(
    Number.of.Genes = n(),
    .groups = "drop"
  )

u = ggplot(plot_data, aes(y = Number.of.Genes, x = region, 
                             group= combo, fill = combo)) +
  geom_bar(stat = 'identity', position = "stack") + theme_bw() + 
  scale_fill_manual(values = c("cadetblue2", "plum2", "lightsalmon", "firebrick2"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
svg("TWAS/Unique_Genes_count.svg", width = 10, height = 6)
u
dev.off()



common.df = plot_data3
unique.df= plot_data[,c(1,2,4)]
colnames(common.df) = c("Region", "Model", "Number.of.Genes")
colnames(unique.df) = c("Region", "Model", "Number.of.Genes")
all.df = rbind(common.df, unique.df)
all.df_noCov = all.df[grep("_Cov", all.df$Model, invert=T),]


a = ggplot(all.df_noCov, aes(y = Number.of.Genes, x = Region, group= Model, fill = Model)) +
  geom_bar(stat = 'identity', position = "dodge") + theme_bw() + 
  scale_fill_manual(values = c("grey", "blue","darkblue", "lightblue", "steelblue",
                               "lightgreen", "darkgreen", "orange", "darkorange", 
                               "brown", "firebrick2"))+
  geom_text(aes(label=Number.of.Genes), size = 4, position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
svg("TWAS/Unique_and_Common_Genes_count.svg", width = 10, height = 6)
a
dev.off()




