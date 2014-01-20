path <- "/Volumes/Governator/ANTISTATELONG/ROIs/Data/dani_scripts"
tbl <- read.table(file.path(path,"linked9to26_20120504_copy.txt"),header=T)

X <- data.frame(
  id = tbl[,which(names(tbl)=="LunaID")],
  age = tbl[,which(names(tbl)=="Age.at.visit")],
  ASpErrCorr = tbl[,which(names(tbl)=="ASpErrCorr")],
	AS.lat.corr.AVG = tbl[,which(names(tbl)=="AS.lat.corr.AVG")],
	AS.lat.errCorr.AVG = tbl[,which(names(tbl)=="AS.lat.errCorr.AVG")],
	VGS.lat.corr.AVG = tbl[,which(names(tbl)=="VGS.lat.corr.AVG")]
)
save(X, file=file.path(path, "analysis20120508", "X"))

Y <- data.frame(numeric(dim(X)[1]))
ynames=c("dlPFC_L", "dlPFC_R", "vlPFC_L", "vlPFC_R", "insula_L", "insula_R", "SEF", "preSMA", "FEF_L", "FEF_R", "putamen_L", "putamen_R", "PPC_L", "PPC_R", "V1_bilat",	"cerebellum_L", "cerebellum_R", "dACC_10")
for(y in 1:length(ynames)) Y[,dim(Y)[2]+1] <- tbl[,which(names(tbl)==ynames[y])]
Y <- Y[,-1]
names(Y) <- ynames
Y <- as.matrix(Y)
save(Y, file=file.path(path, "analysis20120508", "Y"))
