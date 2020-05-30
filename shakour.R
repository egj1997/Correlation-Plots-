#shakour
df <- data.frame(
  Histone = c("H3.1K27M", "H3.3K27M", "Wild-type"),
  Number = c(14, 2, 4)
)
head(df)
library(dplyr)
library(ggsci)
library(wesanderson)
# Compute the position of labels
df <- df %>% 
  arrange(desc(Histone)) %>%
  mutate(prop = Number / sum(df$Number) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

a=ggplot(df, aes(x="", y=prop, fill=Histone)) +
  geom_bar(stat="identity", width=1,color="white") +
  coord_polar("y", start=0)+theme_void()+theme(legend.position="none") +
    geom_text(aes(y = ypos, label = Histone), color = "white", size=6)+scale_fill_simpsons()+ggtitle("G328V/E mutation")
a

df1 <- data.frame(
  Histone = c("H3.1K27M", "H3.3K27M"),
  Number = c(3, 4)
)
df1 <- df1 %>% 
  arrange(desc(Histone)) %>%
  mutate(prop = Number / sum(df$Number) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


b=ggplot(df1, aes(x="", y=prop, fill=Histone)) +
  geom_bar(stat="identity", width=1,color="white") +
  coord_polar("y", start=0)+theme_void()+theme(legend.position="none") +
  geom_text(aes(y = ypos, label = Histone), color = "white", size=6)+scale_fill_locuszoom()+ggtitle("R206H mutation")
b
library(gridExtra)
par(mfrow=c(2,1))
jpeg("Shakour_histone_mutations_piechart.jpeg",res=1200,width = 4, height = 8, units = 'in')

grid.arrange(a,b,nrow=2)

dev.off()
