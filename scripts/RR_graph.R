library(tidyverse); library(data.table)
library(ggsci)

setwd('/home/shanedenecke/Dropbox/Transporter_PK/Cyp6g1_tissues/Jenny_Toxicology')


raw=fread('./27Mar21_Cyp6g1_resistance_ratio_table.csv')
raw$Genotype=gsub('UO','Tubules',raw$Genotype)
raw$Genotype=gsub('Lsp2','Fat Body',raw$Genotype)
raw$Genotype=gsub('MexG','Midgut',raw$Genotype)
raw$Genotype=gsub('HR','All Tissues',raw$Genotype)

raw$Genotype=factor(raw$Genotype,levels=c('Midgut','Tubules','Fat Body','All Tissues'))
raw=arrange(raw,Genotype) %>% mutate(lab=c(rep('A',3),rep('B',3),rep('C',3),rep('D',3)))


gp=ggplot(raw,aes(x=Pesticide,y=ratio,fill=Pesticide))
gp=gp+facet_grid(cols=vars(Genotype))
gp=gp+geom_bar(stat='identity')
gp=gp+geom_errorbar(aes(ymin=lower,ymax=upper))
gp=gp+scale_fill_aaas()
gp=gp+lims(y=c(0,3.5))
gp=gp+labs(y='Resistance Ratio\n',x='\nPesticide')
gp=gp+geom_hline(yintercept = 1,linetype=2)
gp=gp+geom_text(x=.7,y=3.4,size=12,aes(label=lab))
gp=gp+theme_bw()
gp=gp+theme(legend.title=element_text(size=24,face='bold'),
            legend.text=element_text(size=20,face='bold'),
            axis.text.y=element_text(size=20),
            axis.text.x=element_text(size=14,hjust=1,angle=30),
            axis.title=element_text(size=24),
            #panel.grid = element_blank(),
            legend.position='none',
            strip.text = element_text(size=24,hjust = 0.5,face='bold'),
            plot.title = element_text(size=36,hjust = 0.5))
print(gp)

ggsave(plot=gp,filename='./Lc50_RR_figurev2.pdf',device='pdf',width=12,height=6)







gp=ggplot(raw,aes(x=Genotype,y=ratio,fill=Genotype))
gp=gp+facet_grid(cols=vars(Pesticide))
gp=gp+geom_bar(stat='identity')
gp=gp+geom_errorbar(aes(ymin=lower,ymax=upper))
gp=gp+scale_fill_aaas()
gp=gp+labs(y='Resistance Ratio\n',x='Genotype')
gp=gp+geom_hline(yintercept = 1,linetype=2)
gp=gp+theme_bw()
gp=gp+theme(legend.title=element_text(size=24,face='bold'),
            legend.text=element_text(size=20,face='bold'),
            axis.text.y=element_text(size=20),
            axis.text.x=element_text(size=20,hjust=1,angle=30),
            axis.title=element_text(size=24),
            panel.grid = element_blank(),
            legend.position='none',
            strip.text = element_text(size=24,hjust = 0.5),
            plot.title = element_text(size=36,hjust = 0.5))
print(gp)

ggsave(plot=gp,filename='./Lc50_RR_figure.pdf',device='pdf',width=10,height=6)
