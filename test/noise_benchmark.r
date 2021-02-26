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


sa = sample(1:nrow(lexp), 1000)
exp = lexp[sa, ww2[1:500]]

bins = 6


vv = apply(exp,1,sd)

write.table(exp, file="files/testexp.tsv", quote=F, sep="\t")
unlink("files/timeout.tsv")
system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b ",bins," -t 12 -c 2000"))   

rr = read.table("files/timeout.tsv", stringsAsFactors=F, sep="\t")
inter = intersect(colnames(rr), rownames(exp))

exp = exp[inter, ]
rr1 = rr[inter,inter]

cc1 = cor(t(exp[inter,]))


qual = list()

for(i in 1:5){
    mi = c()
    cors = c()
    for(k in 1:1){
        noise = matrix(rnorm(nrow(exp)*ncol(exp)), nrow(exp), ncol(exp))
        for(n in 1:nrow(noise)){
            noise[n,] = noise[n,]*vv[n]*i*1
        }

        exp_noise = noise
        exp_noise[exp_noise < 0] = 0
        cct = cor(t(exp_noise[inter,]))
        write.table(exp_noise, file="files/testexp.tsv", quote=F, sep="\t")

        ctime = Sys.time()
        system(paste0("athena_mi/target/release/athena_mi -i files/testexp.tsv -o files/timeout.tsv -b ",bins," -t 12 -c 2000"))   
        rr = read.table("files/timeout.tsv", stringsAsFactors=F, sep="\t")
        rr = rr[inter, inter]
        qq = cor(unlist(c(rr1)), unlist(c(rr)))
        qqc = cor(unlist(c(cc1)), unlist(c(cct)))
        mi = c(mi, qq)
        cors = c(cors, qqc)
        print(paste0("-----> ",qq," Time elapsed: ", difftime(Sys.time(), ctime, units = "secs")))
    }
    qual[[length(qual)+1]] = cors
    qual[[length(qual)+1]] = mi
}

pdf("figures/noise_benchmark.pdf", 6,6)
boxplot(qual, xlab="", ylab="", cex=1.6, cex.lab=2, cex.axis=2, las=2, border=c("#1a8901", "#ba3412"), col=c("#1a8901","#fc7b5a"), lwd=2)
mtext("threads", side=1, line=5, cex=2.2)
mtext("seconds", side=2, line=5, cex=2.2)
dev.off()


