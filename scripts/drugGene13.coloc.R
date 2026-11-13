#install.packages("dplyr")
#install.packages("data.table")
#install.packages("coloc")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/ieugwasr",force = TRUE)
#devtools::install_github("mrcieu/gwasglue", force = TRUE)
#devtools::install_github("mrcieu/gwasvcf", force = TRUE)


#install.packages("ggplot2")
#install.packages("patchwork")

#引用包
library(dplyr)
library(data.table)
library(coloc)
library(VariantAnnotation)
library(gwasglue)
library(ggplot2)
library(patchwork)
library(geni.plots)

#输入文件
PPH4filter=0.75                          #PPH4打分过滤条件
outcomeFile="finngen_R12_M13_SYSTSLCE.gz"      #结局数据
finnAnnFile="R11_manifest.xls"           #芬兰数据库注释文件
drugFile="drug.IEUid.txt"                #药物数据库ID的文件

#设置工作目录


#读取药物数据库ID的文件
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F)
row.names(drugRT)=drugRT[,"id"]

#获取目录下所有vcf.gz结尾的文件
allVcf=list.files(pattern="*.vcf.gz$")

#对vcf文件循环
outTab=data.frame()
for(eqtlFile in allVcf){
	possibleError=tryCatch({
	#获取基因位置信息
	eqtlID=gsub(".vcf.gz", "", eqtlFile)
	geneChr=drugRT[eqtlID, "Chr"]
	geneStart=drugRT[eqtlID, "Start"]
	geneEnd=drugRT[eqtlID, "End"]
	geneName=drugRT[eqtlID, "Symbol"]

	#读取基因输入文件
	vcfRT <- readVcf(eqtlFile)
	data1=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")

	#挑选暴露数据
	data2 <- data1 %>% dplyr::select("SNP","chr.exposure","pos.exposure","effect_allele.exposure",
	                                 "other_allele.exposure", "eaf.exposure","beta.exposure",
	                                 "se.exposure", "pval.exposure", "samplesize.exposure")
	colnames(data2) <- c("SNP","chrom","Pos","effect_allele","other_allele",
	                     "MAF","beta","se","P","samplesize")
	#整理共定位所需数据
	data3 <- as.data.frame(data2)
	data3$varbeta <- data3$se^2
	data3$z = data3$beta/data3$se
	data3 <- subset(data3, !duplicated(SNP))
	geneData <- data3 %>% filter(chrom==geneChr, Pos>geneStart-100000, Pos<geneEnd+100000) %>% na.omit()
	lead <- geneData %>% dplyr::arrange(P)
	leadPos <- lead$Pos[1]
	QTLdata <-  geneData %>% filter(Pos>leadPos-50000,Pos<leadPos+50000) %>% na.omit()

	#读取结局数据文件
	data0=data.table::fread(outcomeFile, header=T, sep="\t", check.names=F)
	finn_info <- fread(finnAnnFile, data.table=F)
	trait_row <- finn_info[grepl(outcomeFile, finn_info$url),]
	data0$ncase.outcome <- trait_row$num_cases
	data0$ncontrol.outcome <- trait_row$num_controls
	data0$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

	#挑选结局数据
	data1 <- data0 %>% dplyr::select("rsids","#chrom","pos","alt","ref",
	                                 "af_alt","beta","sebeta","pval","samplesize.outcome","ncase.outcome")
	colnames(data1) <- c('SNP','chrom',"pos",'effect_allele','other_allele',
	                     "eaf","beta","se","P","samplesize","number_cases")
	data2 <- as.data.frame(data1)
	data2$varbeta <- data2$se^2
	data2$MAF <- ifelse(data2$eaf<0.5,data2$eaf,1-data2$eaf)
	data3 <- subset(data2, !duplicated(SNP))
	data3$s <- data3$number_cases/data3$samplesize
	data3$z = data3$beta/data3$se
	GWASdata <- data3 %>% na.omit()

	#共定位分析
	sameSNP <- intersect(QTLdata$SNP, GWASdata$SNP)
	QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
	GWASdata <- GWASdata[GWASdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
	coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,
	                                 N=QTLdata$samplesize,MAF=QTLdata$MAF,z = QTLdata$z,
	                                 pvalues=QTLdata$P,type="quant"),
	                   dataset2=list(snp=GWASdata$SNP,beta=GWASdata$beta,varbeta=GWASdata$varbeta,
	                                 N=GWASdata$samplesize,MAF=GWASdata$MAF,z = GWASdata$z,
	                                 pvalues=GWASdata$P,type="cc"))
	result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
	#保存共定位分析的结果
	result$results %>% dplyr::arrange(desc(SNP.PP.H4))
	if(result$summary[6]>=PPH4filter){
		outTab=rbind(outTab, data.frame(ID=eqtlID,
		                                Symbol=geneName,
		                                PP.H0=result$summary[2],
		                                PP.H1=result$summary[3],
		                                PP.H2=result$summary[4],
		                                PP.H3=result$summary[5],
		                                PP.H4=result$summary[6]))

		#===== 用 geni.plots::fig_region_stack 画图 =====
		#1) 合并两个 trait 的 p 值与位置
		plotDat <- merge(
			GWASdata[, c("SNP","pos","P")] %>% dplyr::rename(pos_gwas = pos, P_gwas = P),
			QTLdata[, c("SNP","Pos","P")] %>% dplyr::rename(pos_qtl = Pos, P_qtl = P),
			by = "SNP"
		)
		plotDat <- as.data.frame(plotDat)
		plotDat$pos_gwas <- as.numeric(plotDat$pos_gwas)
		plotDat$P_gwas   <- as.numeric(plotDat$P_gwas)
		plotDat$P_qtl    <- as.numeric(plotDat$P_qtl)

		#2) lead SNP = SNP.PP.H4 最高
		topRow <- result$results %>% dplyr::arrange(desc(SNP.PP.H4)) %>% head(1)
		topSNP <- topRow$snp

		#3) 构造 geni.plots 期望的 assoc 数据框（pos 必须是 integer）
		assoc_df <- data.frame(
			marker   = plotDat$SNP,
			chr      = as.integer(geneChr),
			pos      = as.integer(plotDat$pos_gwas),
			pvalue_1 = plotDat$P_gwas,
			pvalue_2 = plotDat$P_qtl
		)

		#4) 用本地 plink 算的 LD（ld_pos.rds，GRCh37 位置）最近邻匹配（容差 200bp）
		corr <- NULL
		if(file.exists("ld_pos.rds")) {
			ld_pos <- readRDS("ld_pos.rds")
			ld_poss <- as.integer(rownames(ld_pos))
			our_pos_int <- as.integer(assoc_df$pos)
			nearest <- sapply(our_pos_int, function(p) ld_poss[which.min(abs(ld_poss - p))])
			dist <- abs(nearest - our_pos_int)
			keep <- dist <= 200
			cat("DBG local LD fuzzy match:", sum(keep), "/", length(our_pos_int), "max dist:", max(dist[keep]), "\n")
			if(sum(keep) > 5) {
				keep_idx <- which(keep)
				dup <- duplicated(nearest[keep_idx])
				if(any(dup)) {
					cat("DBG dropping", sum(dup), "duplicate LD pos mappings\n")
					keep_idx <- keep_idx[!dup]
				}
				map_pos <- as.character(nearest[keep_idx])
				if(topSNP %in% assoc_df$marker[keep_idx]) {
					corr <- ld_pos[map_pos, map_pos, drop=FALSE]
					rownames(corr) <- colnames(corr) <- assoc_df$marker[keep_idx]
					assoc_df <- assoc_df[keep_idx, ]
				}
			}
		}
		cat("DBG n_assoc:", nrow(assoc_df), "top:", topSNP, "corr:", !is.null(corr), "\n")

		#5) 画图：data 直接传数据框，corr 单独传
		p <- fig_region_stack(
			data        = assoc_df,
			corr        = corr,
			traits      = c("GWAS:finngen_R12_M13_SYSTSLCE", "eQTL：eqtl-a-ENSG00000171860.C3AR1"),
			top_marker  = topSNP,
			r2          = TRUE,
			build       = 38,
			genebar     = FALSE,
			girafe      = FALSE,
			legend      = TRUE
		)
		pdf(file=paste0(eqtlID, "." , geneName, "_v2.pdf"), width=9, height=8)
		print(p)
		dev.off()
	}
	},error=function(e) e)
    if(inherits(possibleError, "error")){next}
}
#输出共定位结果的表格文件
write.csv(outTab, file="coloc.result.csv", row.names=F)
