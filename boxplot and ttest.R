?t.test()
df = read.csv("figure8.csv")
df = df[-26,]
df = df[-26,]
x = df[1:24,3]
y = df[25:35,3]
df$group = "NULL"
df[1:24,4] = "Non-malignant"
df[25:35,4] = "High grade serous carcinoma"
df = cbind(x,y)
t.test(x,y)
boxplot(x,y, xlab = c("x, y"))
?boxplot

library(ggplot2)
ggplot(df, aes(x=group, y=Value, fill=group)) + geom_boxplot() +
  theme(legend.position = "bottom")
       
df = read.csv("figure11.csv")
df$group = "NULL"
df[1:53,4] = "Ovarian Tumor"
df[54:63,4] = "Normal Ovary"
x = df[1:53,3]
y = df[54:63,3]
t.test(x,y)
boxplot(x,y, xlab = c("x, y"))
?boxplot

library(ggplot2)
ggplot(df, aes(x=group, y=Value, fill=group)) + geom_boxplot() +
  theme(legend.position = "bottom")


df = read.csv("figure5_vplot.csv")
df$group = "NULL"
df[1:53,4] = "Ovarian Tumor"
df[54:63,4] = "Normal Ovary"
x = df[1:53,3]
y = df[54:63,3]
t.test(x,y)
boxplot(x,y, xlab = c("x, y"))
?boxplot

df = data
x = data[1:12,3]
y = data[13:69,3]
t.test(x,y)
