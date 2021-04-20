library(tidyverse)
library(data.table)
library(gridExtra)
library(ggsignif)

setwd('/home/shanedenecke/Dropbox/Transporter_PK/Cyp6g1_tissues/Luciferase_GAL4/')



raw=fread('./Final_Luciferase_analysis.csv') %>% 
  #filter((Set=='Fat Body' & Repeat==1) | (Repeat==2 & Set!='Fat Body'))
  #filter(Repeat==2)
  filter(Repeat==1)
raw$`Normalized million RLU`=raw$`Normalized RLU`/1000000
raw$Set=factor(raw$Set,levels=c('Midgut','Tubules','Fat Body'))


gp=ggplot(raw ,aes(x=Cyp6g1,y=`Normalized million RLU`,fill=Cyp6g1))
gp=gp+facet_grid(rows=vars(Set),scales = 'free_y')
gp=gp+lims(y=c(0,25))
gp=gp+stat_summary(fun=mean,position=position_dodge(width=0.95),geom="bar")
gp=gp+geom_dotplot(binaxis = "y",stackdir = "center", colour = "black",position = position_dodge(width = 0.5), dotsize = 1.5)
gp=gp+geom_signif(comparisons=list(c('Hikone-R','Tissue Specific')),test='t.test',y_position=20,tip_length=.2)
gp=gp+labs(x='Genotype',y='Millions RLUs / microgram Protein')
gp=gp+stat_summary(geom="errorbar",colour='black',width=.4,fun.raw="mean_se",size=1)
#gp=gp+scale_y_continuous(expand=expansion(mult = c(0,.5)))
gp=gp+scale_fill_manual(values=c('grey30','red'))
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
ggsave(plot=gp,device='pdf',filename='./Final_luciferase_plot_set1.pdf',width=5,height=10)

l=list()
for(i in unique(raw$Set)){
  sub=raw[Set==i]
  l[[i]]=t.test(sub[Cyp6g1=='Hikone-R']$`Normalized million RLU`,sub[Cyp6g1=='Tissue Specific']$`Normalized million RLU`)$p.value
}


raw %>% group_by(Genotype) %>% summarize(m=mean(`Normalized million RLU`))



gp=ggplot(raw ,aes(x=Set,y=`Normalized million RLU`,fill=Cyp6g1,group=interaction(Set,Cyp6g1)))
#gp=gp+facet_grid(rows=vars(Set),scales = 'free_y')
gp=gp+lims(y=c(0,25))
gp=gp+stat_summary(fun=mean,position=position_dodge(width=0.95),geom="bar")
gp=gp+geom_dotplot(binaxis = "y",stackdir = "center", colour = "black",position = position_dodge(width = 0.95), dotsize = 1.5)
gp=gp+geom_signif(comparisons=list(c('Hikone-R','Tissue Specific')),test='t.test',y_position=20,tip_length=.2)
gp=gp+labs(x='Genotype',y='Millions RLUs / microgram Protein')
gp=gp+stat_summary(geom="errorbar",colour='black',width=.4,position=position_dodge(width=0.95),fun.raw="mean_se",size=1)
#gp=gp+scale_y_continuous(expand=expansion(mult = c(0,.5)))
gp=gp+scale_fill_manual(values=c('grey30','red'))
gp=gp+theme_bw()
gp=gp+theme(legend.title=element_text(size=24,face='bold'),
            legend.text=element_text(size=20,face='bold'),
            axis.text.y=element_text(size=20),
            axis.text.x=element_text(size=20,hjust=1,angle=30),
            axis.title=element_text(size=24),
            panel.grid = element_blank(),
            #legend.position='none',
            strip.text = element_text(size=24,hjust = 0.5),
            plot.title = element_text(size=36,hjust = 0.5))
print(gp)
ggsave(plot=gp,device='pdf',filename='./Final_luciferase_plot_set1.pdf',width=5,height=10)



