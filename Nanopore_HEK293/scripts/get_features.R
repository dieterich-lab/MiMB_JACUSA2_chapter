args <- commandArgs(trailingOnly = TRUE)
#Read preprocessed JACUSA2 output
print("Rebuild matrix of features ...")
Exp1 <-
  read.delim(paste0(args[1], "/data_reformat.txt"),
             as.is = T,
             header = F)

colnames(Exp1) = c(
  "ID",
  "contig",
  "position",
  "call2.score",
  "deletion.score",
  "insertion.score",
  "Base",
  "Anchor",
  "Strand"
)
if (dir.exists(dirname(args[2])) == FALSE) {
  dir.create(dirname(args[2]))
}

Call2 <- tapply(Exp1$call2.score, list(Exp1$ID, Exp1$Anchor), sum)
Call2[is.na(Call2)] <- 0
colnames(Call2) <- paste0("Exp1Call2Score_", colnames(Call2))

Deletion <- tapply(Exp1$deletion.score, list(Exp1$ID, Exp1$Anchor), sum)
Deletion[is.na(Deletion)] <- 0
colnames(Deletion) <- paste0("Exp1DeletionScore_", colnames(Deletion))

Insertion <-
  tapply(Exp1$insertion.score, list(Exp1$ID, Exp1$Anchor), sum)
Insertion[is.na(Insertion)] <- 0
colnames(Insertion) <-
  paste0("Exp1InsertionScore_", colnames(Insertion))

BigTable <-
  merge(
    data.frame(ID = rownames(Call2), Call2),
    data.frame(ID = rownames(Deletion), Deletion),
    by.x = 1,
    by.y = 1
  )
BigTable <-
  merge(BigTable,
        data.frame(ID = rownames(Insertion), Insertion),
        by.x = 1,
        by.y = 1)


motif <-
  read.table(paste0(args[1], "/checkMotif_reformat.txt"),
             as.is = T,
             header = F)
motif[, 1] <- gsub("-", "_", motif[, 1])
motif[, 1] <- gsub("\\(\\_\\)", ":-", motif[, 1])
motif[, 1] <- gsub("\\(\\+\\)", ":\\+", motif[, 1])

BigTable <- merge(BigTable, motif, by.x = 1, by.y = 1)


colnames(BigTable)[ncol(BigTable)] <- "Motif"
rownames(BigTable) <- BigTable[, 1]
BigTable <- BigTable[, -1]

BigTable$DRACH <- rep(0, nrow(BigTable))
BigTable$DRACH[grep("[AGT][AG]AC[ACT]", BigTable$Motif)] <- 1
saveRDS(BigTable, file = args[2])
