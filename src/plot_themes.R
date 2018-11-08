## Graphic (ggplot2) themes
## Script author: David Chen
## Script maintainer: David Chen

require(ggplot2);

myScatterTheme <- theme_classic() + 
  theme(axis.text.x=element_text(size=15,color="black"), axis.title.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=15,color="black"), axis.title.y=element_text(size=15,color="black"),
        title=element_text(size=15,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black")); 
  
mySurvTheme <- theme_classic() +
  theme(axis.text=element_text(size=20,color="black"), 
        axis.title=element_text(size=20,color="black"),
        legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
        legend.box.background=element_rect(colour="black",size=2),
        text=element_text(size=20),
        strip.text.x=element_text(size=15,colour="black",face="bold"));

myBarplotTheme <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=15,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=15, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") ); 

myBoxplotTheme <- theme_classic() +
  theme(axis.text.x=element_text(size=21,color="black"), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=30,color="black"),
        title=element_text(size=15,color="black"));

mySuppBoxplotTheme <- theme_classic() +
  theme(axis.text.x=element_text(face="italic",size=16,color="black"), axis.text.y=element_text(size=16,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=16,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold"));
  
myWaterfallTheme <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=25, color="black"), title=element_text(size=25, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15,color="black"));

