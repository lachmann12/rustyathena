setwd("GitHub/rustyathena")


library("preprocessCore")
library("rhdf5")

gexp = read.table("files/human_10k.tsv", stringsAsFactors=F, sep="\t")
ww = which(rowSums(gexp)/ncol(gexp) > 20)
ww2 = rev(order(colSums(gexp)))
gexp = gexp[ww,]
lexp = normalize.quantiles(log2(data.matrix(gexp)+1))
rownames(lexp) = rownames(gexp)
colnames(lexp) = colnames(gexp)

sa = sample(1:nrow(lexp), 2000)
exp = lexp[sa, ww2[1:1000]]

write.table(exp, file="files/testexp.tsv", quote=F, sep="\t")

ct1 = list()
for(t in 1:6){
    tt = c()
    for(i in 1:10){
        ctime = Sys.time()
        system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b 6 -t ",t," -c 500"))   
        tt = c(tt, difftime(Sys.time(), ctime, units = "secs"))
        print(paste0("----->  Time elapsed: ", difftime(Sys.time(), ctime, units = "secs")))
    }
    ct1[[length(ct1)+1]] = tt
}

ct2 = list()
for(b in 3:12){
    tt = c()
    for(i in 1:10){
        ctime = Sys.time()
        system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b ",b," -t 6 -c 500"))   
        tt = c(tt, difftime(Sys.time(), ctime, units = "secs"))
    }
    ct2[[length(ct2)+1]] = tt
}

names(ct2) = 3:12


ct3 = list()
for(b in 1:10){
    tt = c()
    for(i in 1:10){
        ctime = Sys.time()
        system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b 6 -t 6 -c ", 100*b))   
        tt = c(tt, difftime(Sys.time(), ctime, units = "secs"))
    }
    ct3[[length(ct3)+1]] = tt
}

names(ct3) = (1:10)*100

pdf("figures/speed_benchmark.pdf", 18,6)
#par(mar=c(7,4,1,1))

par(fig=c(0.04,0.37,0.07,1), new=TRUE)
boxplot(ct1, xlab="", ylab="", cex=1.6, cex.lab=2, cex.axis=2, las=2, border="#ba3412", col="#fc7b5a", lwd=2)
mtext("threads", side=1, line=5, cex=2.2)
mtext("seconds", side=2, line=5, cex=2.2)

par(fig=c(0.36,0.69,0.07,1), new=TRUE)
boxplot(ct2, xlab="", ylab="", cex=1.6, cex.lab=2, cex.axis=2, las=2, border="#1a8901", col="#7aef5f", lwd=2)
mtext("bins", side=1, line=5, cex=2.2)

par(fig=c(0.68,1,0.07,1), new=TRUE)
boxplot(ct3, xlab="", ylab="", cex=1.6, cex.lab=2, cex.axis=2, las=2, border="#58008e", col="#c46afc", lwd=2)
mtext("samples", side=1, line=5, cex=2.2)
dev.off()

save(ct1, file="files/ct1.rda")
save(ct2, file="files/ct2.rda")
save(ct3, file="files/ct3.rda")

