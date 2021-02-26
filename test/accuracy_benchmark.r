setwd("GitHub/rustyathena")


library("preprocessCore")
library("rhdf5")

gexp = read.table("files/human_10k.tsv", stringsAsFactors=F, sep="\t")
ww = which(rowMeans(gexp) > 100)
ww2 = rev(order(colSums(gexp)))
gexp = gexp[ww,]
lexp = normalize.quantiles(log2(data.matrix(gexp)+1))
rownames(lexp) = rownames(gexp)
colnames(lexp) = colnames(gexp)


sa = sample(1:nrow(lexp), 4000)
exp = lexp[sa, ww2[1:10000]]

bins = 6


vv = apply(exp,1,sd)

write.table(exp, file="files/testexp.tsv", quote=F, sep="\t")
unlink("files/timeout.tsv")


system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b ",20," -t 12"))   
rr = read.table("files/timeout.tsv", stringsAsFactors=F, sep="\t")

inter = intersect(colnames(rr), rownames(exp))

exp = exp[inter, ]
rr1 = rr[inter,inter]

cc1 = cor(t(exp[inter,]), method="spearman")


cors = list()

for(i in 3:10){
    ctime = Sys.time()
    system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b ",i," -t 12"))   
    rr = read.table("files/timeout.tsv", stringsAsFactors=F, sep="\t")
    rr = rr[inter, inter]
    #qq = cor(unlist(c(rr[upper.tri(rr)])), abs(unlist(c(cc1[upper.tri(cc1)]))))
    qq = diag(cor(rr, abs(cc1)))

    cors[[length(cors)+1]] = qq
    print(paste0("-----> ",mean(qq)," Time elapsed: ", difftime(Sys.time(), ctime, units = "secs")))
}



pdf("figures/accuracy_benchmark.pdf", 7,7)
boxplot(cors)
dev.off()

diag(rr) = NA
diag(cc1) = NA

pdf("figures/test.pdf")
for(i in 1:100){
    em = mean(exp[i,])
    sx = scale(cc1[i,])
    sy = scale(unlist(rr[i,]))
    plot(sx, sy, pch=20, xlab="mutual information", ylab="correlation", main=paste(rownames(rr)[i], format(em, digits=2)))
}
dev.off()

jupyter notebook --port 8000 --no-browser --ip='*' --NotebookApp.token='' --NotebookApp.password='bio'

pdf("figures/accuracy_benchmark.pdf", 7,7)
par(mar=c(6,6,4,4))
boxplot(2:16, cors, xlab="", ylab="", lty=1, type="b", cex=1.6, cex.lab=1.2, cex.axis=1.6, las=2, col="#ba3412", lwd=2)
mtext("bins", side=1, line=4, cex=1.8)
mtext("similarity to abs(correlation)", side=2, line=4, cex=1.8)
dev.off()


pdf("figures/accuracy_benchmark.pdf", 7,7)
boxplot(cors)
dev.off()


