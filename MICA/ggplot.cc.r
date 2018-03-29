library(utils)
library(ggplot2)

s_args = 6
args <- commandArgs()
args <- c(args)
dest_file <- args[s_args + 4]
cc = t(read.table(args[s_args], sep="\t"))
ncols = as.numeric(ncol(cc))
cells = t(cc[1, 2:ncols])

x_tsne <- as.numeric(t(cc[2, 2:ncols]))
y_tsne <- as.numeric(t(cc[3, 2:ncols]))
c <- as.numeric(t(cc[4, 2:ncols])) + 1
uc <- unique(c)
n = as.numeric(length(uc))
df <- data.frame(x_tsne, y_tsne, c)
x_pose = rep(c(0), n)
y_pose = rep(c(0), n)
label = c(1:n)
l_size = rep(c(0), n)
for(i in 1:n){
  x_pose[i] <- sum(x_tsne[c==i]) / length(c[c==i]);
  y_pose[i] <- sum(y_tsne[c==i]) / length(c[c==i]);
  l_size[i] <- paste(c(i, "(", length(x_tsne[c==i]), ")"), collapse = '')
  label[i] <- i;
}
df1 <- data.frame(x_pose, y_pose, label)
p <- ggplot() + 
  geom_point(data=df, aes(x=x_tsne, y=y_tsne, color=as.factor(c)), size=as.numeric(args[s_args+1])) + 
  xlab("MICA-1") + ylab("MICA-2") + 
  guides(col=guide_legend(override.aes = list(size=10), title=paste(c("Clusters\n(", toString(ncols-1), ")"), collapse = '')[1]), ncol=ceiling(n/10)) + 
  geom_text(data=df1, aes(x=x_pose, y=y_pose, label=label), size=as.numeric(args[s_args+2])+4) + 
  scale_color_discrete(labels=l_size) +
  ggtitle(paste(c(args[s_args+3], " (", toString(n), ")"), collapse = '')[1]) + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=10)), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14,face="bold"))
ggsave(dest_file, plot = p, width = 12, height = 12)


