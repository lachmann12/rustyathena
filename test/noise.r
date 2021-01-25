
library("preprocessCore")
library("rhdf5")
h5ls("../Downloads/human_matrix_v9.h5")


gexp = read.table("GitHub/rustyathena/files/human_10k.tsv", stringsAsFactors=F, sep="\t")
ww = which(rowSums(gexp)/ncol(gexp) > 20)
ww2 = rev(order(colSums(gexp)))
gexp = gexp[ww,]
lexp = normalize.quantiles(log2(as.matrix(gexp)+1))
rownames(lexp) = rownames(gexp)
colnames(lexp) = colnames(gexp)

sa = sample(1:nrow(lexp), 5000)
exp = lexp[sa, ww2[1:3000]]

cc = cor(t(exp))

write.table(exp, file="GitHub/rustyathena/files/testexp.tsv", quote=F, sep="\t")

)
ct = c()
for(t in 1:2){
    ctime = Sys.time(
    system(paste0("GitHub/rustyathena/athena_mi/target/release/athena_mi.exe -i GitHub/rustyathena/files/testexp.tsv -o GitHub/rustyathena/files/timeout.tsv -b 12 -t ",," -c 10000"))   
    ct = c(ct, Sys.time()-ctime)
}

system("GitHub/rustyathena/athena_mi/target/release/athena_mi.exe -i GitHub/rustyathena/files/testexp.tsv -o GitHub/rustyathena/files/testout6.tsv -b 6 -t 12 -c 10000")
system("GitHub/rustyathena/athena_mi/target/release/athena_mi.exe -i GitHub/rustyathena/files/testexp.tsv -o GitHub/rustyathena/files/testout20.tsv -b 20 -t 12 -c 10000")


ami = read.table("GitHub/rustyathena/files/testout6.tsv", stringsAsFactors=F, sep="\t", check.names=FALSE)
ami = ami[,rownames(ami)]

ami20 = read.table("GitHub/rustyathena/files/testout20.tsv", stringsAsFactors=F, sep="\t", check.names=FALSE)
ami20 = ami20[,rownames(ami20)]


co = diag(cor(t(cc), t(ami)))


mm = c()
pdf("GitHub/rustyathena/files/mi_cor.pdf", 12, 6)
par(mfrow=c(1,2))
for(k in 1:mm){
    ww = which(cc[k,] < 0.1 & cc[k,] > -0.1 & ami20[k,] > 0.35)
    if(length(ww) > 0) {
        print("wow")
        temp1 = cc[k,]
        temp2 = ami20[k,]
        temp1[k] = NA
        temp2[k] = NA
        plot(temp1, temp2, pch=".", cex=2, main=paste(rownames(cc)[k], k), xlab="Pearson Correlation", ylab="Mutual Information")
        abline(h=0.35, lty=2)
        abline(v=0.1, lty=2)
        abline(v=-0.1, lty=2)
        #mm = c(mm, k)
        ww = which(cc[k,] < 0.1 & cc[k,] > -0.1 & ami20[k,] > 0.35)
        #plot(log(1+unlist(gexp[k,])), log(1+unlist(gexp[ww[1],])), pch=".", cex=2, xlab=rownames(lexp)[k], ylab=rownames(lexp)[ww[1]])
        plot(exp[k,], exp[ww[1],], pch=".", cex=2, xlab=rownames(exp)[k], ylab=rownames(exp)[ww[1]])
    }
}
dev.off()



oo = order(temp1)
tt = temp1[oo]
ot = unlist(temp2[oo])
rsd = rollapply(ot, 100, sd)
rmean = rollapply(ot, 100, mean)

plot(temp1, temp2, pch=".", cex=2, main=paste(rownames(cc)[k], k), xlab="Pearson Correlation", ylab="Mutual Information")
lines(tt[1:length(rmean)], rmean)

ww = which(cc[k,] < 0.1 & cc[k,] > -0.1 & ami[k,] > 0.22)

plot(lexp[k,], lexp[ww[1],], pch=".", cex=3)

ee = colSums(exp)


gexp = read.table("GitHub/rustyathena/files/human_10k.tsv", stringsAsFactors=F, sep="\t")
gexp = gexp[,1:2000]

lexp = normalize.quantiles(log2(as.matrix(gexp)))
coco = cor(t(lexp))




