library(tidyverse)
library(data.table)
library(ggplot2)

setwd('/home/shanedenecke/Dropbox/Transporter_PK/Cyp6g1_tissues/wiggle_summary/')



################### Raw combined #########################

raw=fread('/home/shanedenecke/Dropbox/Transporter_PK/Cyp6g1_tissues/wiggle_summary/all_raw_combined.csv')
times=c('0min','5min','10min','15min','30min','45min','60min','90min','120min')
raw$Time=factor(raw$Time,levels=times)
raw$Set=gsub('(^.+) >.+$','\\1',raw$Genotype)
raw$Expression=gsub('^.+> (.+$)','\\1',raw$Genotype)
raw$individual=paste0(raw$Rep,raw$Location)

flat=spread(raw,key=Time,value=Wiggle_Index)
for(t in times){
  flat[[paste0(t,'__RMR')]]=flat[[t]]/flat[['0min']]
}
rmr.raw=select(flat,Date:individual,contains('RMR'))
rmr.melt=melt(rmr.raw,
         id.vars=c('Date','Genotype','Pesticide','Dose','Rep','Location','Set','individual','Expression'),
         variable.name = 'Time',value.name = 'RMR') %>% 
  mutate(Time=gsub('__RMR','',Time)) %>%
  mutate(Time=factor(Time,levels=times))
           


################## Outliers
ind=select(rmr.melt,Genotype,Pesticide,Dose,Time) %>% unique.data.frame

out.list=list()
for(i in 1:nrow(ind)){
  sub=rmr.melt %>% filter(Genotype==ind[i]$Genotype) %>% filter(Pesticide==ind[i]$Pesticide) %>% filter(Dose==ind[i]$Dose) %>% filter(Time==ind[i]$Time)
  
  outs=boxplot.stats(sub$RMR)$out
  if(length(outs)>0){
    out.list[[paste(ind[i]$Genotype,ind[i]$Pesticide,ind[i]$Dose,ind[i]$Time)]]=sub %>% filter(RMR==outs)
    }
}

out.table=rbindlist(out.list) %>% group_by(Date,Genotype,Pesticide,Dose,individual) %>% summarize(coun=n()) %>% 
  arrange(desc(coun)) %>% filter(coun>=3) %>% select(-coun)
  
  
# rmr.clean=rmr.melt %>% filter(((Date %in% out.table$Date) & 
#                    (Genotype %in% out.table$Genotype) & 
#                    (Pesticide %in% out.table$Pesticide) & 
#                    (Dose %in% out.table$Dose) & 
#                    (individual %in% out.table$individual)))

outliers.final=merge(out.table,rmr.melt,by=colnames(out.table),all.x=T)

rmr.clean=anti_join(rmr.melt,outliers.final)



se=function(x) sd(x)/sqrt(length(x))
rmr.sum=rmr.clean %>% 
  group_by(Date,Genotype,Pesticide,Dose,Time,Expression,Set) %>% 
  summarize(RMR_mean=mean(RMR),RMR_stderr=se(RMR)) %>% data.table()





plot.rmr.sub=rmr.sum %>% filter(
  (Set=='TB' & Dose=='40ppm' & Pesticide=='Bendiocarb') |
  #(Set=='TB' & Dose=='30ppm' & Pesticide=='Imidacloprid') |
  (Set=='TB' & Pesticide=='Imidacloprid' & Date=='06Mar21') |
  (Set=='MG' & Dose=='40ppm' & Pesticide=='Bendiocarb') |
  (Set=='MG' & Dose=='40ppm' & Pesticide=='Imidacloprid') |
  (Set=='FB' & Dose=='40ppm' & Pesticide=='Bendiocarb') |
  (Set=='FB' & Dose=='40ppm' & Pesticide=='Imidacloprid') | 
  (Set=='HR' & Pesticide!='Malathion')
)



plot.rmr.sub$Set=gsub('TB','Tubules',plot.rmr.sub$Set)
plot.rmr.sub$Set=gsub('HR','All Tissues',plot.rmr.sub$Set)
plot.rmr.sub$Set=gsub('MG','Midgut',plot.rmr.sub$Set)
plot.rmr.sub$Set=gsub('FB','Fat Body',plot.rmr.sub$Set)
plot.rmr.sub$Set=factor(plot.rmr.sub$Set,levels=c('Midgut','Tubules','Fat Body','All Tissues'))
plot.rmr.sub$Expression=factor(plot.rmr.sub$Expression,levels=c('NULL','Cyp6g1'))
plot.rmr.sub$Time=as.numeric(gsub('min','',plot.rmr.sub$Time))
                    
v=c()
for(i in 1:nrow(plot.rmr.sub)){
  sub=plot.rmr.sub[i]
  if(sub$Pesticide=="Bendiocarb" & sub$Set=="Midgut"){v[i]="A"}
  if(sub$Pesticide=="Bendiocarb" & sub$Set=="Tubules"){v[i]="B"}
  if(sub$Pesticide=="Bendiocarb" & sub$Set=="Fat Body"){v[i]="C"}
  if(sub$Pesticide=="Bendiocarb" & sub$Set=="All Tissues"){v[i]="D"}
  if(sub$Pesticide=="Imidacloprid" & sub$Set=="Midgut"){v[i]="E"}
  if(sub$Pesticide=="Imidacloprid" & sub$Set=="Tubules"){v[i]="F"}
  if(sub$Pesticide=="Imidacloprid" & sub$Set=="Fat Body"){v[i]="G"}
  if(sub$Pesticide=="Imidacloprid" & sub$Set=="All Tissues"){v[i]="H"}
}
  
plot.rmr.sub$lab=v



gp=ggplot(data=plot.rmr.sub,aes(x=Time,y=RMR_mean,color=Expression))
gp=gp+facet_grid(rows=vars(Pesticide),cols=vars(Set))
gp=gp+scale_y_continuous(breaks=c(0,.25,.5,.75,1,1.25))
gp=gp+scale_color_manual(values=c('black','red'))
gp=gp+labs(x='\nTime (Minutes)',y='Relative Movement Ratio\n')
gp=gp+geom_line(size=1.2)
gp=gp+geom_point(size=4)
gp=gp+geom_errorbar(aes(ymin=RMR_mean-RMR_stderr,ymax=RMR_mean+RMR_stderr))
gp=gp+geom_text(x=1,y=1.5,size=12,aes(label=lab))
gp=gp+theme_bw()
gp=gp+theme(strip.text = element_text(size=24,face='bold'),legend.title = element_text(size=24,face='bold'),
            legend.text = element_text(size=18),
            axis.title = element_text(size=24,face='bold'),
            axis.text=element_text(size=14))
  #panel.grid = element_blank())
print(gp)

ggsave(filename='./Final_Wiggle_Graph.pdf',plot=gp,device='pdf',width=15,height=10)



################################ MALATHION S2



mal.sum=rmr.sum %>% filter(Pesticide=='Malathion' & Set=='HR')
mal.sum$Set=gsub('TB','Tubules',mal.sum$Set)
mal.sum$Set=gsub('HR','All Tissues',mal.sum$Set)
mal.sum$Set=gsub('MG','Midgut',mal.sum$Set)
mal.sum$Set=gsub('FB','Fat Body',mal.sum$Set)
mal.sum$Set=factor(mal.sum$Set,levels=c('Midgut','Tubules','Fat Body','All Tissues'))
mal.sum$Expression=factor(mal.sum$Expression,levels=c('NULL','Cyp6g1'))
mal.sum$Time=as.numeric(gsub('min','',mal.sum$Time))



gp=ggplot(data=mal.sum,aes(x=Time,y=RMR_mean,color=Expression))
gp=gp+facet_grid(rows=vars(Pesticide),cols=vars(Set))
gp=gp+scale_y_continuous(breaks=c(0,.25,.5,.75,1,1.25))
gp=gp+scale_color_manual(values=c('black','red'))
gp=gp+labs(x='Time',y='Relative Movement Ratio')
gp=gp+geom_line(size=1.2)
gp=gp+geom_point(size=4)
gp=gp+geom_errorbar(aes(ymin=RMR_mean-RMR_stderr,ymax=RMR_mean+RMR_stderr))
gp=gp+theme_bw()
gp=gp+theme(strip.text = element_text(size=24,face='bold'),legend.title = element_text(size=24,face='bold'),
            legend.text = element_text(size=18),
            axis.title = element_text(size=18,face='bold'),
            axis.text=element_text(size=14))
#panel.grid = element_blank())
print(gp)


ggsave(filename='./Wiggle_Malathion.pdf',plot=gp,device='pdf',width=8,height=5)






















### #manual filters

plot.rmr=rmr.sum %>% filter(!(Date=='17Mar21' & Set=='MG' & Pesticide=='Imidacloprid')) %>% filter(Pesticide!="Malathion") %>% 
  filter(!(Set=='TB' & Date=='02Mar21'))
plot.rmr$Time=as.numeric(gsub('min','',plot.rmr$Time))


for(i in unique(plot.rmr$Set)){
  sub=plot.rmr %>% filter(Set==i)
  
  
  gp=ggplot(data=sub,aes(x=Time,y=RMR_mean,color=Genotype))
  gp=gp+facet_grid(cols=vars(Pesticide),rows=vars(Dose))
  gp=gp+scale_y_continuous(breaks=c(0,.25,.5,.75,1,1.25))
  gp=gp+ylim(c(0,1.4))
  gp=gp+geom_line()
  gp=gp+geom_point()
  gp=gp+geom_errorbar(aes(ymin=RMR_mean-RMR_stderr,ymax=RMR_mean+RMR_stderr))
  gp=gp+theme_bw()
  print(gp)
  
}




