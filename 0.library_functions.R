# ====================================== #
# Analysis functions used for univariate
# analyses in the Intercell II work
# jamesc@icr.ac.uk, 11th March 2014
# ====================================== #


# This is used at various points where we want to
# know which values in a vector are ≤ -2...
get_ltneg2_per_target <- function(x){
	length(which(x <= -2))
}


compare_deps_across_histologies <- function(
	results,
	marker,
	tissue1,
	tissue2,
	pvalue=0.05,
	rnai_muts,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols,
	file_prefix="compare_deps_across_histologies"
	){
	
	# find the rows for marker in tissue1
	results_to_plot <- which(
		results$marker == marker &
		results$tissue == tissue1 &
		results$PermutationP <= pvalue
		)
	
	filename1 <- paste(file_prefix, tissue1, marker, "depdencies_in", tissue1, ".pdf", sep="_")
	filename2 <- paste(file_prefix, tissue1, marker, "depdencies_in", tissue2, ".pdf", sep="_")
	
	# get targets
	targets <- results[results_to_plot,"target"]
	
	tissue1_kinome_rows <- which(rnai_muts$tissues[,tissue1] == 1)
	tissue2_kinome_rows <- which(rnai_muts$tissues[,tissue2] == 1)
	
	# plot tissue*gene hits in same tissue
	make_box_dot_plots4_col_by_tissue(
		results=as.data.frame(
			results[results_to_plot,]
			),
		zscores=rnai_muts$rnai[tissue1_kinome_rows,],
		mutation.classes=rnai_muts$mut_classes[tissue1_kinome_rows,],
		mutations=rnai_muts$func_muts[tissue1_kinome_rows,],
		exclusions=rnai_muts$all_muts[tissue1_kinome_rows,],
		tissues=rnai_muts$tissues[tissue1_kinome_rows,],
		filename=filename1,
		tissue_pretty_names=tissue_pretty_names,
		tissue_actual_names=tissue_actual_names,
		tissue_cols=tissue_cols
		)
	
	# plot tissue*gene hits in other tissue
	make_box_dot_plots4_col_by_tissue(
		results=as.data.frame(
			results[results_to_plot,]
			),
		zscores=rnai_muts$rnai[tissue2_kinome_rows,],
		mutation.classes=rnai_muts$mut_classes[tissue2_kinome_rows,],
		mutations=rnai_muts$func_muts[tissue2_kinome_rows,],
		exclusions=rnai_muts$all_muts[tissue2_kinome_rows,],
		tissues=rnai_muts$tissues[tissue2_kinome_rows,],
		filename=filename2,
		tissue_pretty_names=tissue_pretty_names,
		tissue_actual_names=tissue_actual_names,
		tissue_cols=tissue_cols
		)
		
}



plot_mutation_barchart <- function(
	mut_classes,
	max_genes=20,
	include_exprn=FALSE,
	pdf_file="mutations_barplot.pdf"
	){
# func start

muts_per_gene_mutd <- NULL
i <- NULL
for(i in 1:ncol(mut_classes)){
	muts_per_gene_mutd <- rbind(
		muts_per_gene_mutd,
		c(
			length(which(mut_classes[,i] == 1)),
			length(which(mut_classes[,i] == 2)),
#			length(which(mut_classes[,i] == 3)),
			length(which(mut_classes[,i] == 4)),
			length(which(mut_classes[,i] == 5))#,
#			length(which(mut_classes[,i] == 6)),
#			length(which(mut_classes[,i] == 7))
			)
		)
}
rownames(muts_per_gene_mutd) <- colnames(mut_classes)
colnames(muts_per_gene_mutd) <- c(
		"missense",
		"truncating",
#		"truncating (hom)",
		"amplified",
		"deleted"#,
#		"overexpressed",
#		"underexpressed"
		)

muts_per_gene_mutd_stripped_names <- NULL
for(i in 1:nrow(muts_per_gene_mutd)){
	muts_per_gene_mutd_stripped_names[i] <- strsplit(rownames(muts_per_gene_mutd)[i], "_", fixed=TRUE)[[1]][1]
}

# need to calculate a sort order for the genes
# based on observed counts.
total_muts_per_gene <- apply(muts_per_gene_mutd,1,sum)
sorted_genes_by_total_muts <- sort(
	total_muts_per_gene,
	index.return=TRUE,
	decreasing=TRUE
	)

# This is the plot of 63 genes with five or more mutants
#pdf(file=pdf_file, width=14, height=7)
pdf(file=pdf_file, width=7, height=5)
par(oma=c(5,0,0,0), mar=c(3,4.2,1,1))

barplot(
	t(muts_per_gene_mutd[sorted_genes_by_total_muts$ix[1:max_genes],]),
	legend.text=c(
		"recurrent\nmissense",
		"truncating",
#		"truncating (hom)",
		"amplification",
		"deletion"#,
#		"overexpressed",
#		"underexpressed"
		),
	col=c("#EE7EA6", "#FF4400", "#C9DD03", "#AF1930"),
#	col=c("#000088", "#FF4400", "#008800", "#AF1930"),
#	col=c("#EE7EA6","#FFD602","#F9A100","#C9DD03","#A71930", "#726E20", "#003D4C"), # pink, yellow, orange, green, red,
# "#EE7EA6" pink
# "#C9DD03" green
# "#A71930" red
# "#F9A100" orange
# "#FFD602" yellow
# "#726E20" olive
# "#003D4C" blue
	space=0,
	las=2,
#	names.arg=rep(NA,30),
	names.arg=muts_per_gene_mutd_stripped_names[sorted_genes_by_total_muts$ix[1:max_genes]],
	ylab="count of mutated cell lines",
	cex.axis=1.5,
	cex.lab=1.5,
	cex.names=1.5
	)
#par(mar=c(4,4.3,0,2))
dev.off()

# func end	
}


read_rnai_drugs_mutations <- function(
	rnai_file,
	drugs_file,
	func_muts_file,
	all_muts_file,
	mut_classes_file,
	tissues_file
	){
	
	rnai <- read.table(
		file=rnai_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	rnai_qn <- t(normalize.quantiles(t(rnai)))
	rownames(rnai_qn) <- rownames(rnai)
	colnames(rnai_qn) <- colnames(rnai)
	
	drugs <- read.table(
	file=drugs_file,
	header=TRUE,
	sep="\t",
	row.names=1
	)
	
	func_muts <- read.table(
		file=func_muts_file,
		sep="\t",
		header=TRUE,
		row.names=1
		)

	all_muts <- read.table(
		file=all_muts_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)

	mut_classes <- read.table(
		file=mut_classes_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	tissues <- read.table(
		file=tissues_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	common_celllines <- intersect(
		rownames(rnai),
		rownames(func_muts)
		)

	common_celllines <- intersect(
		common_celllines,
		rownames(drugs)
		)

	
	i <- NULL
	row.index <- NULL
	rnai_muts_cmn <- NULL
	rnai_qn_muts_cmn <- NULL
	drugs_cmn <- NULL
	func_muts_rnai_cmn <- NULL
	all_muts_rnai_cmn <- NULL
	mut_classes_rnai_cmn <- NULL
	tissues_rnai_cmn <- NULL
	for(i in seq(1:length(common_celllines))){
		# rnai subset
		row.index <- NULL
		row.index <- which(rownames(rnai) == common_celllines[i])
		rnai_muts_cmn <- rbind(
			rnai_muts_cmn,
			rnai[row.index,]
		)
		# rnai_qn subset
		row.index <- NULL
		row.index <- which(rownames(rnai_qn) == common_celllines[i])
		rnai_qn_muts_cmn <- rbind(
			rnai_qn_muts_cmn,
			rnai_qn[row.index,]
		)
		
		# drugs subset
		row.index <- NULL
		row.index <- which(rownames(drugs) == common_celllines[i])
		drugs_cmn <- rbind(
			drugs_cmn,
			drugs[row.index,]
		)
		
		# func_muts subset
		row.index <- NULL
		row.index <- which(rownames(func_muts) == common_celllines[i])
		func_muts_rnai_cmn <- rbind(
			func_muts_rnai_cmn,
			func_muts[row.index,]
		)
		# all_muts subset
		row.index <- NULL
		row.index <- which(rownames(all_muts) == common_celllines[i])
		all_muts_rnai_cmn <- rbind(
			all_muts_rnai_cmn,
			all_muts[row.index,]
		)
		# mut_classes subset
		row.index <- NULL
		row.index <- which(rownames(mut_classes) == common_celllines[i])
		mut_classes_rnai_cmn <- rbind(
			mut_classes_rnai_cmn,
			mut_classes[row.index,]
		)
		# tissues_rnai subset
		row.index <- NULL
		row.index <- which(rownames(tissues) == common_celllines[i])
		tissues_rnai_cmn <- rbind(
			tissues_rnai_cmn,
			tissues[row.index,]
		)
	}
	rownames(rnai_muts_cmn) <- common_celllines
	rownames(func_muts_rnai_cmn) <- common_celllines
	rownames(all_muts_rnai_cmn) <- common_celllines
	rownames(mut_classes_rnai_cmn) <- common_celllines
	rownames(tissues_rnai_cmn) <- common_celllines
	
	# calculate inter-quartile range (iqr) lower fence values
	# as a threshold for sensitivity. Do for the rnai and rnai_qn
	# data sets
	
	rnai_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai)){
		rnai_iqr_stats <- quantile(rnai[,i], na.rm=TRUE)
		rnai_iqr_thresholds[i] <- rnai_iqr_stats[2] - 	((rnai_iqr_stats[4] - rnai_iqr_stats[2]) * 1.5)
	}
	names(rnai_iqr_thresholds) <- colnames(rnai)
	
	rnai_qn_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai_qn)){
		rnai_qn_iqr_stats <- quantile(rnai_qn[,i], na.rm=TRUE)
		rnai_qn_iqr_thresholds[i] <- rnai_qn_iqr_stats[2] - 	((rnai_qn_iqr_stats[4] - rnai_qn_iqr_stats[2]) * 1.5)
	}
	names(rnai_qn_iqr_thresholds) <- colnames(rnai_qn)


	return(
		list(
			rnai=rnai_muts_cmn,
			rnai_qn=rnai_qn_muts_cmn,
			drugs=drugs_cmn,
			func_muts=func_muts_rnai_cmn,
			all_muts=all_muts_rnai_cmn,
			mut_classes=mut_classes_rnai_cmn,
			rnai_qn_iqr_thresholds=rnai_qn_iqr_thresholds,
			rnai_iqr_thresholds=rnai_iqr_thresholds,
			tissues=tissues_rnai_cmn
			)
		)	
}





read_drugs_mutations <- function(
	drugs_file,
	func_muts_file,
	all_muts_file,
	mut_classes_file,
	mutation_data_tissues_file,
	tissues_file
	){
	
	
	drugs <- read.table(
		file=drugs_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	func_muts <- read.table(
		file=func_muts_file,
		sep="\t",
		header=TRUE,
		row.names=1
		)

	all_muts <- read.table(
		file=all_muts_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)

	mut_classes <- read.table(
		file=mut_classes_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)

	mutation_data_tissues <- read.table(
		file=tissues_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	common_celllines <- intersect(
		rownames(drugs),
		rownames(func_muts)
		)

	
	i <- NULL
	row.index <- NULL
	drugs_cmn <- NULL
	func_muts_cmn <- NULL
	all_muts_cmn <- NULL
	mut_classes_cmn <- NULL
	tissues_cmn <- NULL
	
	for(i in seq(1:length(common_celllines))){
		
		# drugs subset
		row.index <- NULL
		row.index <- which(rownames(drugs) == common_celllines[i])
		drugs_cmn <- rbind(
			drugs_cmn,
			drugs[row.index,]
		)
		
		# func_muts subset
		row.index <- NULL
		row.index <- which(rownames(func_muts) == common_celllines[i])
		func_muts_cmn <- rbind(
			func_muts_cmn,
			func_muts[row.index,]
		)
		
		# all_muts subset
		row.index <- NULL
		row.index <- which(rownames(all_muts) == common_celllines[i])
		all_muts_cmn <- rbind(
			all_muts_cmn,
			all_muts[row.index,]
		)
		
		# mut_classes subset
		row.index <- NULL
		row.index <- which(rownames(mut_classes) == common_celllines[i])
		mut_classes_cmn <- rbind(
			mut_classes_cmn,
			mut_classes[row.index,]
		)
		
		# tissues subset
		row.index <- NULL
		row.index <- which(rownames(mutation_data_tissues) == common_celllines[i])
		tissues_cmn <- rbind(
			tissues_cmn,
			mutation_data_tissues[row.index,]
		)
		
	}
	rownames(drugs_cmn) <- common_celllines
	rownames(func_muts_cmn) <- common_celllines
	rownames(all_muts_cmn) <- common_celllines
	rownames(mut_classes_cmn) <- common_celllines
	rownames(tissues_cmn) <- common_celllines
	
	return(
		list(
			drugs=drugs_cmn,
			func_muts=func_muts_cmn,
			all_muts=all_muts_cmn,
			mut_classes=mut_classes_cmn,
			tissues=tissues_cmn
			)
		)	
}

#####


make_mini_box_dot_plots4_col_by_tissue <- function(
	results,
	zscores,
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols
	){

	pdf(file=filename, width=2.5, height=3)

	par(bty="n", tcl=-0.2, mai=c(0.75, 0.7, 0.1, 0.1)) # turn off boxes for plots
	# adding ylbias=-0.5 is not helpful
	
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]
			
			# make a factor with three levels:
			# 		wt,
			#		non-recurrent mutant
			# 		recurrent mutant
			# use for boxplot and  stripchart x-axis
			
			# start by setting all cell lines to wt
			wt_mut_grps_strings <- rep(
				"wt",
				times=length(mutations[,results$marker[i]])
				)
			
			# set the non-recurrent mutations
			wt_mut_grps_strings[which(exclusions[,results$marker[i]] == 1)] <- "non-rec. mut."
			
			# set the recurrent/functional mutations
			wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "rec. mut."
			
			wt_grp_rows <- which(wt_mut_grps_strings == "wt")
			nonfunc_mut_grp_rows <- which(wt_mut_grps_strings == "non-rec. mut.")
			func_mut_grp_rows <- which(wt_mut_grps_strings == "rec. mut.")
			
						
		# boxplot based on all data (wt and mut groups)
			boxplot(
				zscores[wt_grp_rows,results$target[i]],
				zscores[func_mut_grp_rows,results$target[i]],
				pch="",
#				sub=paste(marker_gene, "status"),
#				ylab=paste(target_gene, "z-score"),
				names=c("wt", "mutant"),
				cex.axis=1.5
				)
			abline(
				-2,0,col="red",lty=2
				)
			
			mtext(paste(marker_gene, "status"), 1, line=2, cex=1.5)
			mtext(paste(target_gene, "Z-score"), 2, line=2.2, cex=1.5)

#			print(summary(wt_mut_grps))
			
			# points for each tissue type
			j <- NULL
			for(j in 1:length(tissue_actual_names)){
				tissue <- tissue_actual_names[j]				
				wt_rows_by_tissue <- which(
					wt_mut_grps_strings == "wt" &
					tissues[,tissue] == 1
					)
				mutant_rows_by_tissue <- which(
					wt_mut_grps_strings == "rec. mut." &
					tissues[,tissue] == 1
					)
				
				
				if(length(wt_rows_by_tissue) > 0){
					# plot at 1
					points(
						jitter(rep(1,times=length(wt_rows_by_tissue)), amount=0.33),
						zscores[wt_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19,
						cex=1.5
						)
				}
				if(length(mutant_rows_by_tissue) > 0){
					# plot at 2
					points(
						jitter(rep(2,times=length(mutant_rows_by_tissue)), amount=0.33),
						zscores[mutant_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19,
						cex=1.5
						)
				}
			}		
		}
	}
	
	dev.off()

}





#####


make_box_dot_plots4_col_by_tissue <- function(
	results,
	zscores,
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols,
	response_type="Z-score"
	){

	pdf(file=filename, width=2, height=3)

	par(bty="n", tcl=-0.2, mai=c(0.75, 0.7, 0.1, 0.1)) # turn off boxes for plots
	# adding ylbias=-0.5 is not helpful
	
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]
			
			# make a factor with three levels:
			# 		wt,
			#		non-recurrent mutant
			# 		recurrent mutant
			# use for boxplot and  stripchart x-axis
			
			# start by setting all cell lines to wt
			wt_mut_grps_strings <- rep(
				"wt",
				times=length(mutations[,results$marker[i]])
				)
			
			# set the non-recurrent mutations
			wt_mut_grps_strings[which(exclusions[,results$marker[i]] == 1)] <- "non-rec. mut."
			
			# set the recurrent/functional mutations
			wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "rec. mut."
			
			wt_grp_rows <- which(wt_mut_grps_strings == "wt")
			nonfunc_mut_grp_rows <- which(wt_mut_grps_strings == "non-rec. mut.")
			func_mut_grp_rows <- which(wt_mut_grps_strings == "rec. mut.")
			
						
		# boxplot based on all data (wt and mut groups)
			boxplot(
				zscores[wt_grp_rows,results$target[i]],
				zscores[func_mut_grp_rows,results$target[i]],
				pch="",
#				sub=paste(marker_gene, "status"),
#				ylab=paste(target_gene, "z-score"),
				names=c("wt", "mutant")
				)
			abline(
				-2,0,col="red",lty=2
				)
			
			mtext(paste(marker_gene, "status"), 1, line=2)
			mtext(paste(target_gene, response_type), 2, line=2)

#			print(summary(wt_mut_grps))
			
			# points for each tissue type
			j <- NULL
			for(j in 1:length(tissue_actual_names)){
				tissue <- tissue_actual_names[j]				
				wt_rows_by_tissue <- which(
					wt_mut_grps_strings == "wt" &
					tissues[,tissue] == 1
					)
				mutant_rows_by_tissue <- which(
					wt_mut_grps_strings == "rec. mut." &
					tissues[,tissue] == 1
					)
				
				
				if(length(wt_rows_by_tissue) > 0){
					# plot at 1
					points(
						jitter(rep(1,times=length(wt_rows_by_tissue)), amount=0.33),
						zscores[wt_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19
						)
				}
				if(length(mutant_rows_by_tissue) > 0){
					# plot at 2
					points(
						jitter(rep(2,times=length(mutant_rows_by_tissue)), amount=0.33),
						zscores[mutant_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19
						)
				}
			}		
		}
	}
	
	dev.off()

}





make_box_dot_plots3_col_by_tissue <- function(
	results,
	zscores,
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols
	){

	pdf(file=filename, width=4, height=4)

	par(bty="n", tcl=-0.2) # turn off boxes for plots

	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]
			
			# make a factor with three levels:
			# 		wt,
			#		non-recurrent mutant
			# 		recurrent mutant
			# use for boxplot and  stripchart x-axis
			
			# start by setting all cell lines to wt
			wt_mut_grps_strings <- rep(
				"wt",
				times=length(mutations[,results$marker[i]])
				)
			
			# set the non-recurrent mutations
			wt_mut_grps_strings[which(exclusions[,results$marker[i]] == 1)] <- "non-rec. mut."
			
			# set the recurrent/functional mutations
			wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "rec. mut."
			
			wt_grp_rows <- which(wt_mut_grps_strings == "wt")
			nonfunc_mut_grp_rows <- which(wt_mut_grps_strings == "non-rec. mut.")
			func_mut_grp_rows <- which(wt_mut_grps_strings == "rec. mut.")
			
						
		# boxplot based on all data (wt and mut groups)
			boxplot(
				zscores[wt_grp_rows,results$target[i]],
				zscores[func_mut_grp_rows,results$target[i]],
				pch="",
				names=c("wt", "mutant"),
				sub=paste(marker_gene, "status"),
				ylab=paste(target_gene, "z-score")
				)
			abline(
				-2,0,col="red",lty=2
				)

#			print(summary(wt_mut_grps))
			
			# points for each tissue type
			j <- NULL
			for(j in 1:length(tissue_actual_names)){
				tissue <- tissue_actual_names[j]				
				wt_rows_by_tissue <- which(
					wt_mut_grps_strings == "wt" &
					tissues[,tissue] == 1
					)
				mutant_rows_by_tissue <- which(
					wt_mut_grps_strings == "rec. mut." &
					tissues[,tissue] == 1
					)
				
				
				if(length(wt_rows_by_tissue) > 0){
					# plot at 1
					points(
						jitter(rep(1,times=length(wt_rows_by_tissue)), amount=0.33),
						zscores[wt_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19
						)
				}
				if(length(mutant_rows_by_tissue) > 0){
					# plot at 2
					points(
						jitter(rep(2,times=length(mutant_rows_by_tissue)), amount=0.33),
						zscores[mutant_rows_by_tissue,results$target[i]],
						col=legend_col[j],
						pch=19
						)
				}
			}		
		}
	}
	
	dev.off()

}


#
# make a function very similar to make_box_plots3
# to handle expression data ~ mutation status
# require explicit definition of marker and target
#


make_exprn_box_dot_plots3_col_by_tissue <- function(
	marker,
	target,
#	results,
	zscores, # this is expression
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols,
	tissue
	){

	pdf(file=filename, width=2, height=3)

	par(bty="n", tcl=-0.2, mai=c(0.75, 0.7, 0.1, 0.1)) # turn off boxes for plots

	marker_gene <- strsplit(marker, "_")[[1]][1]
	target_gene <- strsplit(target, "_")[[1]][1]
	
	# start by setting all cell lines to wt
	wt_mut_grps_strings <- rep(
		"wt",
		times=length(mutations[,marker])
		)
	
	# set the non-recurrent mutations
	wt_mut_grps_strings[which(exclusions[,marker] == 1)] <- "non-rec. mut."
	
	# set the recurrent/functional mutations
	wt_mut_grps_strings[which(mutations[,marker] == 1)] <- "rec. mut."

	wt_grp_rows <- which(wt_mut_grps_strings == "wt")
	nonfunc_mut_grp_rows <- which(wt_mut_grps_strings == "non-rec. mut.")
	func_mut_grp_rows <- which(wt_mut_grps_strings == "rec. mut.")
				
   # boxplot based on all data (wt and mut groups)
 	boxplot(
		zscores[wt_grp_rows,target],
		zscores[func_mut_grp_rows,target],
		pch="",
#		sub=paste(marker_gene, "status"),
#		ylab=paste(target_gene, "expression"),
		names=c("wt", "mutant")
		)
	abline(
		-2,0,col="red",lty=2
		)

		mtext(paste(marker_gene, "status"), 1, line=2)
		mtext(paste(target_gene, "expression"), 2, line=2)


	
	# points for each tissue type
#	j <- NULL
#	for(j in 1:length(tissue_actual_names)){
#		tissue <- tissue_actual_names[j]		# tissue is now set in the functions arguments
				
		wt_rows_by_tissue <- which(
			wt_mut_grps_strings == "wt" &
			tissues[,tissue] == 1
			)
		mutant_rows_by_tissue <- which(
			wt_mut_grps_strings == "rec. mut." &
			tissues[,tissue] == 1
			)
		
		if(length(wt_rows_by_tissue) > 0){
			# plot at 1
			points(
				jitter(rep(1,times=length(wt_rows_by_tissue)), amount=0.33),
				zscores[wt_rows_by_tissue,target],
				col=tissue_cols,
				pch=19
				)
		}
		if(length(mutant_rows_by_tissue) > 0){
			# plot at 2
			points(
				jitter(rep(2,times=length(mutant_rows_by_tissue)), amount=0.33),
				zscores[mutant_rows_by_tissue,target],
				col=tissue_cols,
				pch=19
				)
		}
#	}		

	dev.off()
}








make_box_dot_plots2_col_by_tissue <- function(
	results,
	zscores,
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols
	){

	pdf(file=filename, width=6, height=4)
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]
			
			# make a factor with three levels:
			# 		wt,
			#		non-recurrent mutant
			# 		recurrent mutant
			# use for boxplot and  stripchart x-axis
			
			# start by setting all cell lines to wt
			wt_mut_grps_strings <- rep(
				"wt",
				times=length(mutations[,results$marker[i]])
				)
			
			# set the non-recurrent mutations
			wt_mut_grps_strings[which(exclusions[,results$marker[i]] == 1)] <- "non-rec. mut."
			
			# set the recurrent/functional mutations
			wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "rec. mut."
			
			# create a factor from the string.
			wt_mut_grps <- factor(wt_mut_grps_strings, levels=c("wt", "non-rec. mut.", "rec. mut."))
			
		# boxplot based on all data
			boxplot(
				zscores[,results$target[i]] ~ wt_mut_grps,
				#zscores[,results$target[i]] ~ mutations[,results$marker[i]],
				pch="",
				names=c("wt", "mutant", "recurrent mutant"),
				sub=paste(marker_gene, "status"),
				ylab=paste(target_gene, "z-score")
				)
			abline(
				-2,0,col="red",lty=2
				)

#			print(summary(wt_mut_grps))
			
			# points for each tissue type
			j <- NULL
			for(j in 1:length(tissue_actual_names)){
				tissue <- tissue_actual_names[j]
				rows_to_plot <- which(tissues[,tissue] == 1)
				
				if(length(rows_to_plot) < 1){
					next
				}
				print(tissue)
				print(length(rows_to_plot))
				
				wt_mut_grps_sub <- wt_mut_grps[rows_to_plot]
				zscores_sub <- NULL # if we have only one cell line we need to fix the data matrix to have proper dimensions
				
				# skip if no rows to plot
				if(length(rows_to_plot) < 1){
					next
				}
				
				if(length(rows_to_plot) == 1){
					zscores_sub <- matrix(zscores[rows_to_plot,], nrow=length(rows_to_plot), ncol=ncol(zscores), byrow=FALSE, dimnames=list(rownames(zscores)[rows_to_plot], colnames(zscores)))
				}else{
						zscores_sub <- as.matrix(zscores[rows_to_plot,], nrow=length(rows_to_plot), ncol=ncol(zscores))
				}
				
				mutation.classes_sub <- mutation.classes[rows_to_plot,]
				mutations_sub <- mutations[rows_to_plot,]
				exclusions_sub <- exclusions[rows_to_plot,]
				tissues_sub <- tissues[rows_to_plot,]
				
				stripchart(
					zscores_sub[,results$target[i]] ~ wt_mut_grps_sub,
					pch=19,
					col=legend_col[j],
					vertical=TRUE,
					add=TRUE,
					method="jitter",
					jitter=0.3
					)

			}
			
		}
	
	}
	
	dev.off()

}







run_kruskal_by_tissue <- function(
	kinome_file=kinome_file,
	tissue_factors_file=tissue_factors_file
	){
	
	kinome <- read.table(
		file=kinome_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	tissue_factors <- read.table(
		file=tissue_factors_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	factor_levels <- levels(tissue_factors$tissue)
	kw_results <- NULL
	
	i <- NULL
	for(i in 1:ncol(kinome)){
		temp_kw_result <- kruskal.test(
			kinome[,i]~tissue_factors$tissue
			)
		
		dunn.pairs.to.print <- ""
		if(temp_kw_result$p.value <= 0.05){
			temp_dunn_result <- dunn.test(
				kinome[,i],tissue_factors$tissue
				)
			
			pair_labels <- NULL
			j <- NULL # count BONE to OVARY
			for(j in 2:length(factor_levels)){
				k <- 1 # count BREAST to PANCREAS
				while(k < j){
				#for(k in (j-1):(length(factor_levels)-1)){
					pair_labels <- c(
						pair_labels,
						paste(
							factor_levels[j],
							factor_levels[k],
							sep=":"
							)
						)
					k <- k + 1
				}
			}
			
			dunn.pairs.sig <- pair_labels[which(
				temp_dunn_result$P <= 0.05
				)]
			
			for(pair in dunn.pairs.sig){
				dunn.pairs.to.print <- paste(
					dunn.pairs.to.print,
					pair,
					sep=","
					)
			}
		}
		
		
		kw_results <- rbind(
			kw_results,
			c(
				"tissue.factors",
				colnames(kinome)[i],
				temp_kw_result$parameter,
				temp_kw_result$statistic,
				temp_kw_result$p.value,
				dunn.pairs.to.print
				)
			)


	}
	colnames(kw_results) <- c(
		"marker",
		"target",
		"df",
		"kw.chisq.stat",
		"kw.pvalue",
		"dunn.sig.tissue.pairs"
		)
	return(kw_results)
}





test_kinome_with_self_cnv_by_tissue <- function(x){
	
	tissue_types <- colnames(x$tissues)
	
	uv_results_bytissue <- NULL
	tissue <- NULL
	for(tissue in tissue_types){
		
		cellline_count <- sum(
			x$tissues[,tissue]
			)
		
		if(cellline_count < 5){
			next
		}
		
		tissue_rows <- which(
			x$tissues[,tissue] == 1
			)
		
		temp_results <- NULL
		temp_results <- test_kinome_with_self_cnv(
			zscores=x$rnai[tissue_rows,],
			mutations=x$func_muts[tissue_rows,],
			all_variants=x$all_muts[tissue_rows,],
			sensitivity_thresholds=x$rnai_iqr_thresholds
			)
		
		if(is.null(nrow(temp_results))){
			print(paste("Skipping ", tissue, " - no results", sep=""))
			next
		}
		
		temp_results <- cbind(
			temp_results,
			rep(tissue, times=nrow(temp_results))
			)
		
		uv_results_bytissue <- rbind(
			uv_results_bytissue,
			temp_results
			)
	}
	colnames(
		uv_results_bytissue
		)[ncol(
			uv_results_bytissue
			)
			] <- "tissue"
	return(uv_results_bytissue)
	
	
}


run_univariate_test_bytissue <- function(x){
	
	tissue_types <- colnames(x$tissues)
	
	
	uv_results_bytissue <- NULL
	tissue <- NULL
	for(tissue in tissue_types){
		
		cellline_count <- sum(
			x$tissues[,tissue]
			)
		
		if(cellline_count < 5){
			next
		}
		
		tissue_rows <- which(
			x$tissues[,tissue] == 1
			)
		
		temp_results <- NULL
		temp_results <- run_univariate_tests(
			zscores=x$rnai[tissue_rows,],
			mutations=x$func_muts[tissue_rows,],
			all_variants=x$all_muts[tissue_rows,],
			sensitivity_thresholds=x$rnai_iqr_thresholds
			)
		
		if(is.null(nrow(temp_results))){
			print(paste("Skipping ", tissue, " - no results", sep=""))
			next
		}
		
		temp_results <- cbind(
			temp_results,
			rep(tissue, times=nrow(temp_results))
			)
		
		uv_results_bytissue <- rbind(
			uv_results_bytissue,
			temp_results
			)
	}
	colnames(
		uv_results_bytissue
		)[ncol(
			uv_results_bytissue
			)
			] <- "tissue"
	return(uv_results_bytissue)
}



define_sensitising_targets <- function(
	rnai,
	histogram_file,
	scatter_file
	){
	
	rnai <- as.matrix(rnai)
	
	# calc Pearson 2nd skew coeffs
	# count the number of cell lines below z == -2
	rnai_skews <- NULL
	rnai_num_lt_neg2 <- NULL
	i <- NULL
	for(i in 1:ncol(rnai)){
		rnai_skews[i] <- 3 * (mean(rnai[,i], na.rm=TRUE) - median(rnai[,i], na.rm=TRUE)) / sd(rnai[,i], na.rm=TRUE)
		rnai_num_lt_neg2[i] <- length(which(rnai[,i] <= -2))
	}
	names(rnai_skews) <- colnames(rnai)
	names(rnai_num_lt_neg2) <- colnames(rnai)
	
	# pick the sensitising targets in the kinome
	rnai_sensitising_targets <- NULL
	rnai_sensitising_targets_info <- NULL
	i <- NULL
	for(i in 1:length(rnai_skews)){
		if(rnai_skews[i] <= median(rnai_skews, na.rm=TRUE) && rnai_num_lt_neg2[i] >= 3){
			rnai_sensitising_targets <- c(
				rnai_sensitising_targets,
				names(rnai_skews)[i]
				)
			rnai_sensitising_targets_info <- rbind(
				rnai_sensitising_targets_info,
				c(
					rnai_skews[i],
					rnai_num_lt_neg2[i]
					)
				)
		}
	}
	rownames(rnai_sensitising_targets_info) <- rnai_sensitising_targets
	colnames(rnai_sensitising_targets_info) <- c("skewness", "number_lt_neg2")
	
	
	# plot histograms for selected targets
	pdf(histogram_file)
	i <- NULL
	for(i in 1:length(rnai_sensitising_targets)){
		hist(
			rnai[,rnai_sensitising_targets[i]],
			breaks=50,
			main=rnai_sensitising_targets[i],
			xlab="z-score",
			col="red"
			)
	}
	dev.off()
	
	# Save the set of sensitising targets
	# not done for this function...
	
	# visualise the selected targets
	pdf(file=scatter_file, width=4, height=4)
	plot(
		NULL,
		NULL,
		xlim=c(-1.3,0.5),
		ylim=c(0,120),
		xlab="skewness",
		ylab="number of sensitive cell lines"
		)
	rect(
		-1.3,5,-0.37,120,
		density=10,
		col="red",
		border=NA
		)
	points(
		rnai_skews,
		rnai_num_lt_neg2,
		pch=19,
		col=rgb(0,0,0,0.25),
		)
	dev.off()
	
	return(rnai_skews)
}


# Function to read in data, find the intersecting
# rownames and return a list of dataframes with 
# the common cell lines.
# Added 10th Dec 2014 to stop the data preprocessing
# getting out of hand.
read_rnai_mutations <- function(
	rnai_file,
	func_muts_file,
	all_muts_file,
	mut_classes_file,
	tissues_file
	){

	
	rnai <- read.table(
		file=rnai_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	rnai_qn <- t(normalize.quantiles(t(rnai)))
	rownames(rnai_qn) <- rownames(rnai)
	colnames(rnai_qn) <- colnames(rnai)
	
	func_muts <- read.table(
		file=func_muts_file,
		sep="\t",
		header=TRUE,
		row.names=1
		)

	all_muts <- read.table(
		file=all_muts_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)

	mut_classes <- read.table(
		file=mut_classes_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	tissues <- read.table(
		file=tissues_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	common_celllines <- intersect(
		rownames(rnai),
		rownames(func_muts)
		)

	common_celllines <- intersect(
		common_celllines,
		rownames(tissues)
		)


	
	i <- NULL
	row.index <- NULL
	rnai_muts_cmn <- NULL
	rnai_qn_muts_cmn <- NULL
	func_muts_rnai_cmn <- NULL
	all_muts_rnai_cmn <- NULL
	mut_classes_rnai_cmn <- NULL
	tissues_rnai_cmn <- NULL
	for(i in seq(1:length(common_celllines))){
		# rnai subset
		row.index <- NULL
		row.index <- which(rownames(rnai) == common_celllines[i])
		rnai_muts_cmn <- rbind(
			rnai_muts_cmn,
			rnai[row.index,]
		)
		# rnai_qn subset
		row.index <- NULL
		row.index <- which(rownames(rnai_qn) == common_celllines[i])
		rnai_qn_muts_cmn <- rbind(
			rnai_qn_muts_cmn,
			rnai_qn[row.index,]
		)
		# func_muts subset
		row.index <- NULL
		row.index <- which(rownames(func_muts) == common_celllines[i])
		func_muts_rnai_cmn <- rbind(
			func_muts_rnai_cmn,
			func_muts[row.index,]
		)
		# all_muts subset
		row.index <- NULL
		row.index <- which(rownames(all_muts) == common_celllines[i])
		all_muts_rnai_cmn <- rbind(
			all_muts_rnai_cmn,
			all_muts[row.index,]
		)
		# mut_classes subset
		row.index <- NULL
		row.index <- which(rownames(mut_classes) == common_celllines[i])
		mut_classes_rnai_cmn <- rbind(
			mut_classes_rnai_cmn,
			mut_classes[row.index,]
		)
		# tissues_rnai subset
		row.index <- NULL
		row.index <- which(rownames(tissues) == common_celllines[i])
		tissues_rnai_cmn <- rbind(
			tissues_rnai_cmn,
			tissues[row.index,]
		)
	}
	rownames(rnai_muts_cmn) <- common_celllines
	rownames(func_muts_rnai_cmn) <- common_celllines
	rownames(all_muts_rnai_cmn) <- common_celllines
	rownames(mut_classes_rnai_cmn) <- common_celllines
	rownames(tissues_rnai_cmn) <- common_celllines
	
	# calculate inter-quartile range (iqr) lower fence values
	# as a threshold for sensitivity. Do for the rnai and rnai_qn
	# data sets
	
	rnai_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai)){
		rnai_iqr_stats <- quantile(rnai[,i], na.rm=TRUE)
		rnai_iqr_thresholds[i] <- rnai_iqr_stats[2] - 	((rnai_iqr_stats[4] - rnai_iqr_stats[2]) * 1.5)
	}
	names(rnai_iqr_thresholds) <- colnames(rnai)
	
	rnai_qn_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai_qn)){
		rnai_qn_iqr_stats <- quantile(rnai_qn[,i], na.rm=TRUE)
		rnai_qn_iqr_thresholds[i] <- rnai_qn_iqr_stats[2] - 	((rnai_qn_iqr_stats[4] - rnai_qn_iqr_stats[2]) * 1.5)
	}
	names(rnai_qn_iqr_thresholds) <- colnames(rnai_qn)


	return(
		list(
			rnai=rnai_muts_cmn,
			rnai_qn=rnai_qn_muts_cmn,
			func_muts=func_muts_rnai_cmn,
			all_muts=all_muts_rnai_cmn,
			mut_classes=mut_classes_rnai_cmn,
			rnai_qn_iqr_thresholds=rnai_qn_iqr_thresholds,
			rnai_iqr_thresholds=rnai_iqr_thresholds,
			tissues=tissues_rnai_cmn
			)
		)
}



#
# stripchart significant target_expression * mutant gene associations
#
plot_tissue_exprn_by_mut_status <- function(
	expression,
	mutations,
	tissues,
	marker,
	target
	){
	
	# thow away markers or targets we don't need
	expression <- expression[,target]
	mutations <- mutations[,marker]
		
	# This is difficult. If you stripchart the mutants
	# and non-mutants separately, the dots are out of
	# alignment because there are a different number of
	# factors in the mutant and non-mutant tissue types.
	# A horrible hack is to add two of every level to the
	# ends of the data where one is set to mutant and the
	# other set of non-mutant. These will need to be given
	# expression values that are off the scale so they
	# don't get plotted :( 
	
	num_tissue_levels <- length(levels(as.factor(tissues)))

	# get the max value of expression before we append
	expression_max <- max(expression)
	expression_min <- min(expression)
	
	expression <- c(
		expression,
		rep(
			1000000,
			times=(num_tissue_levels*2)
			)
		)

	start_index_added_levels <- length(expression) - ((2*num_tissue_levels)-1) # gate post error...
	end_index_added_levels <- length(expression)

	names(
		expression
		)[start_index_added_levels:end_index_added_levels] <- c(
			levels(as.factor(tissues)),levels(as.factor(tissues))
			)

	tissues <- c(
		tissues,
		levels(as.factor(tissues)),
		levels(as.factor(tissues))
		)
	mutations <- c(
		mutations,
		rep(1, times=(num_tissue_levels)),
		rep(0, times=(num_tissue_levels))
		)
	
	stripchart(
		log2(
			expression[which(mutations == 0)]
			)~
		as.factor(
			tissues[which(mutations == 0)]
			),
		vertical=TRUE,
		pch=19,
		col=rgb(0,0,0,0.5),
		method="jitter",
		jitter=0.2,
		las=2,
		ylim=c(log2(expression_min),log2(expression_max)),
		ylab=paste(target, " log2 expression"),
		main=paste("mutated gene: ", marker)
		)

	stripchart(
		log2(
			expression[which(mutations == 1)]
			)~
		as.factor(
			tissues[which(mutations == 1)]
			),
		vertical=TRUE,
		pch=19,
		col=rgb(1,0,0,0.5),
		method="jitter",
		jitter=0.2,
		las=2,
		add=TRUE
		)
}



#
# Note that this function should be given a subset of statistically
# significant results to process. Additional filters within the
# function are made:
#		at least three mutant lines
#
# modified on 11th March to add the option to set the lower limit of the y-axis 
# defaults to 'auto'
make_box_dot_plots2 <- function(results, zscores, mutation.classes, mutations, exclusions, tissues, filename, ymin="auto"){

	pdf(file=filename, width=4, height=4)
	
	
	auto_determine_y <- FALSE
	if(ymin == "auto"){
		auto_determine_y = TRUE
	}
	
	
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]

			tissue <- NULL
			rows_to_plot <- NULL
			if(!is.null(results$tissue[i])){
				tissue <- results$tissue[i]
				rows_to_plot <- which(tissues[,tissue] == 1)
			}else{
				tissue <- "PanCan"
				rows_to_plot <- 1:nrow(zscores)
			}
			
			# need to subset the data sets here based on tissue
			if(tissue == "PanCan"){
				zscores_sub <- zscores
				mutation.classes_sub <- mutation.classes
				mutations_sub <- mutations
				exclusions_sub <- exclusions
				tissues_sub <- tissues
			}else{
				zscores_sub <- zscores[rows_to_plot,]
				mutation.classes_sub <- mutation.classes[rows_to_plot,]
				mutations_sub <- mutations[rows_to_plot,]
				exclusions_sub <- exclusions[rows_to_plot,]
				tissues_sub <- tissues[rows_to_plot,]
			}
			
			
			#print(rows_to_plot)
			
			if(auto_determine_y){
				ymin=min(zscores_sub[,results$target[i]], na.rm=TRUE)
			}
			
			boxplot(
				zscores_sub[,results$target[i]] ~ mutations_sub[,results$marker[i]],
				pch="",
				names=c("wt", "mutant"),
				sub=paste(marker_gene, "status"),
				ylab=paste(target_gene, "z-score"),
				main=paste(tissue, "\n", round(results$mptest.p[i], 3), sep=""),
				ylim=c(ymin,2)
				)
	
			abline(
				-2,0,col="red",lty=2
				)
			
			# plot for amps/dels/missense/truncs/wt as different colours
			# mutation.classes.cmn

			und_exp_celllines <- which(mutation.classes_sub[,results$marker[i]] == 7)
			ovr_exp_celllines <- which(mutation.classes_sub[,results$marker[i]] == 6)
	
			del_celllines <- which(mutation.classes_sub[,results$marker[i]] == 5)
			amp_celllines <- which(mutation.classes_sub[,results$marker[i]] == 4)
			trunc_celllines <- which(mutation.classes_sub[,results$marker[i]] == 3 | mutation.classes_sub[,results$marker[i]] == 2)
			missense_celllines <- which(mutation.classes_sub[,results$marker[i]] == 1)
			wt_celllines <- which(mutation.classes_sub[,results$marker[i]] == 0)


			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[und_exp_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[und_exp_celllines],results$target[i]],
				pch=19,
				col=rgb(0.00,0.24,0.30,0.5) # under expressed
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[ovr_exp_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[ovr_exp_celllines],results$target[i]],
				pch=19,
				col=rgb(0.45,0.43,0.13,0.5) # over expressed
				)
			
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[del_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[del_celllines],results$target[i]],
				pch=19,
				col=rgb(0.65,0.1,0.18,0.5) #dels
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[amp_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[amp_celllines],results$target[i]],
				pch=19,
				col=rgb(0.79,0.86,0.01,0.75) # amps
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[trunc_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[trunc_celllines],results$target[i]],
				pch=19,
				col=rgb(0.97,0.63,0,0.5) # truncs
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[missense_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[missense_celllines],results$target[i]],
				pch=19,
				col=rgb(0.93,0.49,0.65,0.5) # missense
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[wt_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[wt_celllines],results$target[i]],
				pch=19,
				col=rgb(0.38,0.39,0.39,0.5) # wt
				)
	
		}
	
	}
	
	dev.off()

}





#
# Deprecated... Use make_box_dot_plots2
#
make_box_dot_plots <- function(results, zscores, mutation.classes, mutations, exclusions, tissues, filename){

	pdf(file=filename, width=4, height=4)
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

#		if(results$mptest.p[i] <= 0.05 & results$med.grpA.med.grpB[i] <= -0.3 & results$nA[i] >= 4){
#		if(results$mptest.p[i] <= 0.05 & results$med.grpA.med.grpB[i] <= -0.3 & results$nA[i] >= 4){
#		if(results$mptest.p[i] <= 0.05 & results$med.grpA.med.grpB[i] <= -0.3 & results$nA[i] >= 4){
		if(results$mptest.p[i] <= 0.05 & results$spearman.r[i] <= -0.25 & results$nA[i] >= 3 & results$min.grpA[i] <= -2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]

			tissue <- NULL
			rows_to_plot <- NULL
			if(!is.null(results$tissue[i])){
				tissue <- results$tissue[i]
				rows_to_plot <- which(tissues[,tissue] == 1)
			}else{
				tissue <- "PanCan"
				rows_to_plot <- 1:nrow(kinome)
			}
			
			# need to subset the data sets here based on tissue
			if(tissue == "PanCan"){
				zscores_sub <- zscores
				mutation.classes_sub <- mutation.classes
				mutations_sub <- mutations
				exclusions_sub <- exclusions
				tissues_sub <- tissues
			}else{
				zscores_sub <- zscores[rows_to_plot,]
				mutation.classes_sub <- mutation.classes[rows_to_plot,]
				mutations_sub <- mutations[rows_to_plot,]
				exclusions_sub <- exclusions[rows_to_plot,]
				tissues_sub <- tissues[rows_to_plot,]
			}
			
			
			#print(rows_to_plot)
			
			boxplot(
				zscores_sub[,results$target[i]] ~ mutations_sub[,results$marker[i]],
				pch="",
				names=c("wt", "mutant"),
				sub=paste(marker_gene, "status"),
				ylab=paste(target_gene, "z-score"),
				main=tissue
				)
	
			abline(
				-2,0,col="red",lty=2
				)
			
			# plot for amps/dels/missense/truncs/wt as different colours
			# mutation.classes.cmn
	
			del_celllines <- which(mutation.classes_sub[,results$marker[i]] == 5)
			amp_celllines <- which(mutation.classes_sub[,results$marker[i]] == 4)
			trunc_celllines <- which(mutation.classes_sub[,results$marker[i]] == 3 | mutations_sub[,results$marker[i]] == 2)
			missense_celllines <- which(mutation.classes_sub[,results$marker[i]] == 1)
			wt_celllines <- which(mutation.classes_sub[,results$marker[i]] == 0)
			
			
			
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[del_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[del_celllines],results$target[i]],
				pch=19,
				col=rgb(0.65,0.1,0.18,0.5) #dels
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[amp_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[amp_celllines],results$target[i]],
				pch=19,
				col=rgb(0.79,0.86,0.01,0.75) # amps
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[trunc_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[trunc_celllines],results$target[i]],
				pch=19,
				col=rgb(0.97,0.63,0,0.5) # truncs
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[missense_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[missense_celllines],results$target[i]],
				pch=19,
				col=rgb(0.93,0.49,0.65,0.5) # missense
				)
			points(
				jitter(mutations_sub[rownames(mutation.classes_sub)[wt_celllines],results$marker[i]], amount=0.2)+1,
				zscores_sub[rownames(mutation.classes_sub)[wt_celllines],results$target[i]],
				pch=19,
				col=rgb(0.38,0.39,0.39,0.5) # wt
				)
	
		}
	
	}
	
	dev.off()

}


classify_mutants_by_threshold <- function(zscores, mutations, all_variants, thresholds){
	
	# find grpA and grpB using mutations and exclusions
	# count number of each group below threshold
	# write out:
	# marker, target, nA, nB, threshold, num_sens_grpA, num_sens_grpB, propn_sens_grpA, propn_sens_grpB
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		
		#gene <- strsplit(
		#	colnames(mutations)[i],
		#	"_"
		#	)[[1]][1]
		gene <- colnames(mutations)[i]
		
		
		# grpB includes cell lines with no reported mutations at all
		# in gene...
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3 as we are never going to
		# consider anything based on n=2
		if(length(grpA) < 3 | length(grpB) < 5){
			next
		}
		
		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){

			threshold <- thresholds[colnames(zscores)[j]]
			num_sens_grpA <- length(which(zscores[grpA,j] <= threshold))
			num_sens_grpB <- length(which(zscores[grpB,j] <= threshold))
			num_res_grpA <- length(which(zscores[grpA,j] > threshold))
			num_res_grpB <- length(which(zscores[grpB,j] > threshold))
			
			if(num_sens_grpA == 0){
				next
			}
#			if(sum(num_sens_grpA,num_sens_grpB,num_res_grpA,num_res_grpB) == 0){
#				next
#			}
			if(sum(num_sens_grpA,num_sens_grpB) == 0){
				next
			}
			
			cont.table <- as.matrix(rbind(c(num_sens_grpA,num_res_grpA),c(num_sens_grpB,num_res_grpB)))
			colnames(cont.table) <- c("sens", "res")
			rownames(cont.table) <- c("mut", "wt")
			fisher_result <- fisher.test(cont.table)
			
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			sensitivity <- num_sens_grpA / nA 
			specificity <- num_res_grpB / nB
#			propn_sens_grpA <- num_sens_grpA / nA
#			propn_sens_grpB <- num_sens_grpB / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			
		results <- rbind(
				results,
				c(
					marker,
					target,
					nA,
					nB,
					threshold,
					med.grpA,
					med.grpB,
					med.diff,
					num_sens_grpA,
					num_sens_grpB,
#					propn_sens_grpA,
#					propn_sens_grpB,
					sensitivity,
					specificity,
					min.grpA,
					min.grpB,
					fisher_result$p.value,
					fisher_result$estimate
					)
				)
			}
		}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
#		"fraction.grpA.sens",
#		"fraction.grpB.sens",
		"sensitivity",
		"specificity",
		"min.grpA",
		"min.grpB",
		"fisher.exact.p",
		"odds.ratio"
		)	
	
	return(results)

}



#
# Note - this is deprecated by the classify_mutants_by_threshold()
# function which does the same thing with a list of thresholds
# passed as a parameter... 
#
classify_mutants_by_iqd <- function(zscores, mutations, all_variants){
	
	# find grpA and grpB using mutations and exclusions
	# count number of each group below threshold
	# write out:
	# marker, target, nA, nB, threshold, num_sens_grpA, num_sens_grpB, propn_sens_grpA, propn_sens_grpB
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		
		#gene <- strsplit(
		#	colnames(mutations)[i],
		#	"_"
		#	)[[1]][1]
		gene <- colnames(mutations)[i]
		
		
		# grpB includes cell lines with no reported mutations at all
		# in gene...
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3 as we are never going to
		# consider anything based on n=2
		if(length(grpA) < 3 | length(grpB) < 5){
			next
		}
		
		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){

			quantile_stats <- quantile(zscores[,j], na.rm=TRUE)
			quantile_threshold <-quantile_stats[2] - (quantile_stats[4] - quantile_stats[2])
			num_sens_grpA <- length(which(zscores[grpA,j] <= quantile_threshold))
			num_sens_grpB <- length(which(zscores[grpB,j] <= quantile_threshold))
			num_res_grpA <- length(which(zscores[grpA,j] >= quantile_threshold))
			num_res_grpB <- length(which(zscores[grpB,j] >= quantile_threshold))
			
			if(num_sens_grpA == 0){
				next
			}
			if(sum(num_sens_grpA,num_sens_grpB,num_res_grpA,num_res_grpB) == 0){
				next
			}
			
			cont.table <- as.matrix(rbind(c(num_sens_grpA,num_res_grpA),c(num_sens_grpB,num_res_grpB)))
			colnames(cont.table) <- c("sens", "res")
			rownames(cont.table) <- c("mut", "wt")
			fisher_result <- fisher.test(cont.table)
			
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			propn_sens_grpA <- num_sens_grpA / nA
			propn_sens_grpB <- num_sens_grpB / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			
		results <- rbind(
				results,
				c(
					marker,
					target,
					nA,
					nB,
					quantile_threshold,
					med.grpA,
					med.grpB,
					med.diff,
					num_sens_grpA,
					num_sens_grpB,
					propn_sens_grpA,
					propn_sens_grpB,
					min.grpA,
					min.grpB,
					fisher_result$p.value,
					fisher_result$estimate
					)
				)
			}
		}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
		"fraction.grpA.sens",
		"fraction.grpB.sens",
		"min.grpA",
		"min.grpB",
		"fisher.exact.p",
		"odds.ratio"
		)	
	
	return(results)

}




plot_auc_and_mutations <- function(
	kinome,
	mutations,
	threshold,
	target,
	tree_genes, # list of genes used in tree
	file="auc_and_mutations_plot.pdf",
	bottom_margin=20
	){
	
	# missing values (in z-scores) are a massive pain
	# for this reason, I'm joining the rnai and mutation
	# data and na.omit()ing the lot from the start...
	
	target_mutations <- cbind(
		kinome[,target],
		mutations[,tree_genes]
		)
	target_mutations <- na.omit(target_mutations)
	colnames(target_mutations)[1] <- target
	
	pathway_status <- apply(target_mutations[,-1],1,max)
	target_mutations <- cbind(pathway_status, target_mutations)
	colnames(target_mutations)[1] <- "pathway_status"

	
	# get the sort order for cell lines by z-score
	zscore_sort_order <- sort(
		target_mutations[,target],
		index.return=TRUE
		)
	
	# any cell lines with missing values for z-scores
	# above will have been removed from zscore_sort_order
	# need to be able to exclude these cell lines from the 
	# mutation data too.
	
	
	celllines_ordered <- rownames(target_mutations)[zscore_sort_order$ix]
	
	celllines_ordered_axis_labels <- sub("_.+", "", celllines_ordered, perl=TRUE)
	
	# z-score heatmap colour palette
	breaks=seq(0.5, 1, by=0.01) 
	breaks=append(breaks,0,0)
	mycol <- colorpanel(n=length(breaks)-1,low="cyan",high="black")
	
	
	pdf(
		file=file,
		width=14,
		height=5
		)
	
	layout(matrix(c(1,2,3,3,3,3,3,3,3,3),10,1,byrow=TRUE), respect=FALSE)
	par(oma=c(0,0,0.1,0))
	

	# add z-score colour bar
	par(mar=c(0.5,9,0.5,6))
	
	image(
		as.matrix(target_mutations[celllines_ordered,target]),
		xaxt="n",
		yaxt="n",
		col=mycol,
		breaks=breaks
		)
	axis(
		side=2,
		at=0.5,
		labels=paste(strsplit(target,"_")[[1]][1], "\nAUC", sep=""),
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)

	# add pathway status
	par(mar=c(0.5,9,0.5,6))
	
	image(
		as.matrix(target_mutations[celllines_ordered,"pathway_status"]),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#303030")
		)
	axis(
		side=2,
		at=0.5,
		labels="Pathway\nstatus",
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	
	# decide bottom margin size based on number of marker genes
	#bottom_margin_size <- round(30 - (length(tree_genes_to_print) * 1.429), digits=0)
	
	
	# add mutation status by gene
	par(mar=c(bottom_margin,9,0.5,6))
	image(
		as.matrix(
			target_mutations[
				celllines_ordered,
				tree_genes
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A7A7A7")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(target_mutations)-1)),
		labels=celllines_ordered_axis_labels,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	# clean up tree_gene names for printing
	tree_genes_to_print <- NULL
	i <- NULL
	for(i in 1:length(tree_genes)){
		tree_genes_to_print[i] <- strsplit(tree_genes[i],"_")[[1]][1]
	}
	
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(tree_genes)-1)),
	#	at=seq(from=0, to=1, by=1/(14-1)),
		labels=tree_genes_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
	
}






plot_z_and_mutations3 <- function(
	kinome,
	mutations,
	threshold,
	target,
	tree_genes, # list of genes used in tree
	file="z_and_mutations_plot.pdf",
	bottom_margin=20,
	right_margin=6
	){
	
	# missing values (in z-scores) are a massive pain
	# for this reason, I'm joining the rnai and mutation
	# data and na.omit()ing the lot from the start...
	
	target_mutations <- cbind(
		kinome[,target],
		mutations[,tree_genes]
		)
	target_mutations <- na.omit(target_mutations)
	colnames(target_mutations)[1] <- target
	
	pathway_status <- apply(target_mutations[,-1],1,max)
	target_mutations <- cbind(pathway_status, target_mutations)
	colnames(target_mutations)[1] <- "pathway_status"

	
	# get the sort order for cell lines by z-score
	zscore_sort_order <- sort(
		target_mutations[,target],
		index.return=TRUE
		)
	
	# any cell lines with missing values for z-scores
	# above will have been removed from zscore_sort_order
	# need to be able to exclude these cell lines from the 
	# mutation data too.
	
	
	celllines_ordered <- rownames(target_mutations)[zscore_sort_order$ix]
	
	celllines_ordered_axis_labels <- sub("_.+", "", celllines_ordered, perl=TRUE)
	
	# z-score heatmap colour palette
	breaks=seq(-4, 4, by=0.2) 
	breaks=append(breaks, 100)
	breaks=append(breaks, -100, 0)
	mycol <- colorpanel(n=length(breaks)-1,low="cyan",mid="black",high="yellow")
	
	
	pdf(
		file=file,
		width=14,
		height=5
		)
	
	layout(matrix(c(1,2,3,3,3,3,3,3,3,3),10,1,byrow=TRUE), respect=FALSE)
	par(oma=c(0,0,0.1,0))
	

	# add z-score colour bar
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,target]),
		xaxt="n",
		yaxt="n",
		col=mycol,
		breaks=breaks
		)
	axis(
		side=2,
		at=0.5,
		labels=paste("si", strsplit(target,"_")[[1]][1], "\nz-score", sep=""),
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)

	# add pathway status
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,"pathway_status"]),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#303030")
		)
	axis(
		side=2,
		at=0.5,
		labels="Pathway\nstatus",
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	
	# decide bottom margin size based on number of marker genes
	#bottom_margin_size <- round(30 - (length(tree_genes_to_print) * 1.429), digits=0)
	
	
	# add mutation status by gene
	par(mar=c(bottom_margin,9,0.5,right_margin))
	image(
		as.matrix(
			target_mutations[
				celllines_ordered,
				tree_genes
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A7A7A7")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(target_mutations)-1)),
		labels=celllines_ordered_axis_labels,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	# clean up tree_gene names for printing
	tree_genes_to_print <- NULL
	i <- NULL
	for(i in 1:length(tree_genes)){
		tree_genes_to_print[i] <- strsplit(tree_genes[i],"_")[[1]][1]
	}
	
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(tree_genes)-1)),
	#	at=seq(from=0, to=1, by=1/(14-1)),
		labels=tree_genes_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
	
}







plot_z_and_mutations2 <- function(
	kinome,
	mutations,
	threshold,
	target,
	tree_genes, # list of genes used in tree
	file="z_and_mutations_plot.pdf"
	){
	
	# missing values (in z-scores) are a massive pain
	# for this reason, I'm joining the rnai and mutation
	# data and na.omit()ing the lot from the start...
	
	target_mutations <- cbind(
		kinome[,target],
		mutations
		)
	target_mutations <- na.omit(target_mutations)
	colnames(target_mutations)[1] <- target
	
	# get the sort order for cell lines by z-score
	zscore_sort_order <- sort(
		target_mutations[,target],
		index.return=TRUE
		)
	
	# any cell lines with missing values for z-scores
	# above will have been removed from zscore_sort_order
	# need to be able to exclude these cell lines from the 
	# mutation data too.
	
	
	celllines_ordered <- rownames(target_mutations)[zscore_sort_order$ix]
		
	pdf(
		file=file,
		width=14,
		height=7
		)
	
	layout(matrix(c(1,2,2),3,1,byrow=TRUE), respect=FALSE)
	par(oma=c(0,0,0,0))
	par(mar=c(0.5,6,0,2))
	plot(
		target_mutations[celllines_ordered,target],
		xaxt="n",
		xlab="",
		ylab=paste("si", strsplit(target,"_")[[1]][1], " z-score", sep=""),
		#ylim=c(-10,4),
		pch=19,
		bty="n",
		cex=1.5,
		cex.axis=1.5,
		cex.lab=1.5
		)
	abline(
		threshold,
		0,
		lty=2
		)
	
	par(mar=c(20,9,0.5,6))
	image(
		as.matrix(
			target_mutations[
				celllines_ordered,
				tree_genes
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A71930")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(target_mutations)-1)),
		labels=celllines_ordered,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	# clean up tree_gene names for printing
	tree_genes_to_print <- NULL
	i <- NULL
	for(i in 1:length(tree_genes)){
		tree_genes_to_print[i] <- strsplit(tree_genes[i],"_")[[1]][1]
	}
	
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(tree_genes)-1)),
	#	at=seq(from=0, to=1, by=1/(14-1)),
		labels=tree_genes_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
	
}






plot_z_and_mutations <- function(
	kinome,
	mutations,
	threshold,
	target,
	sensitive_lines,
	tree_genes, # list of genes used in tree
	file="z_and_mutations_plot.pdf"
	){
	
	celllines_ordered <- names(
	sort(kinome[,target], na.last=TRUE)
	)
	
	# mark up each ordered cell line with '*' for
	# predicted sensitive and ' ' for not
	celllines_ordered_to_print <- celllines_ordered
	i <- NULL
	for(i in 1:length(celllines_ordered_to_print)){
		line_is_sig <- celllines_ordered_to_print[i] %in% sensitive_lines
		if(line_is_sig){
			celllines_ordered_to_print[i] <- paste(
				celllines_ordered_to_print[i],
				"*",
				sep=""
				)
		}else{
			celllines_ordered_to_print[i] <- paste(
				celllines_ordered_to_print[i],
				" ",
				sep=""
				)
		}
	}

	
	pdf(
		file=file,
		width=14,
		height=7
		)
	
	layout(matrix(c(1,2,2),3,1,byrow=TRUE), respect=FALSE)
	par(oma=c(0,0,0,0))
	par(mar=c(0.5,6,0,2))
	plot(
		kinome[celllines_ordered,target],
		xaxt="n",
		xlab="",
		ylab=paste("si", target, " z-score", sep=""),
		ylim=c(-10,4),
		pch=19,
		bty="n",
		cex=1.5,
		cex.axis=1.5,
		cex.lab=1.5
		)
	abline(
		threshold,
		0,
		lty=2
		)
	
	par(mar=c(20,9,0.5,6))
	image(
		as.matrix(
			mutations[
				celllines_ordered,
				tree_genes
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A71930")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(mutations)-1)),
		labels=celllines_ordered_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(tree_genes)-1)),
	#	at=seq(from=0, to=1, by=1/(14-1)),
		labels=tree_genes,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
	
}


run_PLS_analysis <- function(
	zscores,
	mutations,
	outfile="PLS_analysis_plots.pdf"
	){
	
	# To hold:
	# min RMSEP
	# num sig comps
	# best predictors
	this.pls.result <- NULL
	
	pdf(file=outfile)
	
	full.pls <- plsr(
			zscores ~ mutations,
			ncomp=10,
			scale=FALSE,
			validation="LOO",
			jackknife=TRUE
			)
	
	
	this.pls.result <- c(
		min(RMSEP(full.pls)$val[1,1,]),
		RMSEP(full.pls)$comps[
			which(
				RMSEP(
				full.pls)$val[1,1,] == min(
					RMSEP(full.pls)$val[1,1,])
					)
				]
		)
	
	plot(
		RMSEP(full.pls),
		legendpos="topright"
		)
	
	fewercomps.pls <- plsr(
			zscores ~ mutations,
			ncomp=2,
			scale=FALSE,
			jackknife=TRUE,
			validation="LOO"
			)
	
	plot(
		fewercomps.pls,
		ncomp=2,
		asp=1,
		line=TRUE
		)
	
	fewercomps.pls.jk <- jack.test(fewercomps.pls, ncomp=2, use.mean=FALSE)
	fewercomps.sig.preds <- names(fewercomps.pls.jk$pvalues[which(fewercomps.pls.jk$pvalues < 0.1),1,1])

# write this into the results table
#	fewercomps.sig.preds

	this.pls.result <- c(
	this.pls.result,
	paste(list(fewercomps.sig.preds), sep="_")
	)
	
	barplot(fewercomps.pls.jk$coefficients[1:30,1,1], las=2)
	
	fewercomps.pls.sig.preds <- plsr(
			zscores ~ mutations[,fewercomps.sig.preds],
			ncomp=2,
			scale=FALSE,
			validation="LOO",
			jackknife=TRUE
			)
	
	plot(
		fewercomps.pls.sig.preds,
		ncomp=2,
		asp=1,
		line=TRUE
		)
		dev.off()
		
		return(this.pls.result)

} # End run_PLS_analysis



run_simple_PLS_analysis <- function(
	zscores,
	mutations
	){
	
	# To hold:
	# min RMSEP
	# num sig comps
	# best predictors
	# min z
	# nominally sensitive cell lines
	this.pls.result <- NULL
		
	full.pls <- plsr(
			zscores ~ mutations,
			ncomp=10,
			scale=FALSE,
			validation="LOO",
			jackknife=TRUE
			)
	
	
#	max_comp_to_use <- min(RMSEP(full.pls)$val[1,1,])
	
	full.pls.jk <- jack.test(full.pls, ncomp=2, use.mean=FALSE)
	full.sig.preds <- names(full.pls.jk$pvalues[which(full.pls.jk$pvalues < 0.1),1,1])

	min_z <- min(zscores)
	
	num_sens <- length(which(zscores == 1))

	this.pls.result <- c(
		min(RMSEP(full.pls)$val[1,1,]),
		RMSEP(full.pls)$comps[
			which(
				RMSEP(
				full.pls)$val[1,1,] == min(
					RMSEP(full.pls)$val[1,1,])
					)
				],
		min_z,
		num_sens,
		paste(list(full.sig.preds), sep="_")
		)
		
		return(this.pls.result)

} # End run_PLS_analysis



run_simple_dichotomised_PLS_analysis <- function(
	zscores,
	mutations
	){
	
	# To hold:
	# min RMSEP
	# num sig comps
	# best predictors
	# min z
	# nominally sensitive cell lines
	this.pls.result <- NULL
		
	full.pls <- plsr(
			zscores ~ mutations,
			ncomp=10,
			scale=FALSE,
			validation="LOO",
			jackknife=TRUE
			)
	
	
#	max_comp_to_use <- min(RMSEP(full.pls)$val[1,1,])
	
	full.pls.jk <- jack.test(
		full.pls,
		ncomp=2,
		use.mean=FALSE
		)
	full.sig.preds <- names(
		full.pls.jk$pvalues[which(
			full.pls.jk$pvalues < 0.1
			),1,1]
		)

#	min_z <- min(zscores)
	
	num_sens <- length(which(zscores == 1))

		
	# This gives the predicted values after 2 comps
	# full.pls$validation$pred[,,2]
	# need to compare this to the actual values
	# to calculate the number of correct preds
	TP <- length(which(zscores == 1 & full.pls$validation$pred[,,2] > 0.5))
	FP <- length(which(zscores == 0 & full.pls$validation$pred[,,2] > 0.5))
	TN <- length(which(zscores == 0 & full.pls$validation$pred[,,2] <= 0.5))
	FN <- length(which(zscores == 1 & full.pls$validation$pred[,,2] <= 0.5))


	this.pls.result <- c(
		min(RMSEP(full.pls)$val[1,1,]), # ie. max_comp_to_use
		#max_comp_to_use,
		RMSEP(full.pls)$comps[
			which(
				RMSEP(
				full.pls)$val[1,1,] == min(
					RMSEP(full.pls)$val[1,1,])
					)
				],
#		min_z,
		num_sens,
		TP,
		FP,
		TN,
		FN,
		paste(list(full.sig.preds), sep="_")
		)

	names(this.pls.result) <- c(
		"min RMSEP",
		"sig comps",
#		"min z",
		"num sens",
		"TP",
		"FP",
		"TN",
		"FN",
		"best predictors"
		)

		return(this.pls.result)

} # End run_PLS_analysis





make_heatmaps_isogenics <- function(
	results,
	zscores,
	mutations,
	all_variants,
	outfile="isogenics_zscore_heatmaps_by_mutated_genes.pdf"
	){	
	pdf(file=outfile)
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		gene <- strsplit(
			colnames(mutations)[i],
			"_"
			)[[1]][1]
		grpB <- which(all_variants[,gene] == 0)
		grpAB <- c(grpA,grpB)
		gene.hits.sens <- results[which(results[,"marker"] == gene & as.numeric(results[,"spearman.p"]) <= 0.05 & as.numeric(results[,"med.grpA.med.grpB"]) < -0), 2]
		gene.hits.res <- results[which(results[,"marker"] == gene & as.numeric(results[,"spearman.p"]) <= 0.05 & as.numeric(results[,"med.grpA.med.grpB"]) >= 0), 2]
	
		# change the test considered below to Spearman from 
		# mptest. mptest can not acheive small p-values with
		# only three samples per group.
		gene.hits.sens.sort <- sort(
			as.numeric(results[which(results[,"marker"] == gene &
			as.numeric(results[,"spearman.p"]) <= 0.05 &
			as.numeric(results[,"med.grpA.med.grpB"]) < -1.5), "med.grpA.med.grpB"]),
			index.return=TRUE, decreasing=TRUE)
	
		gene.hits.res.sort <- sort(
			as.numeric(results[which(results[,"marker"] == gene &
			as.numeric(results[,"spearman.p"]) <= 0.05 &
			as.numeric(results[,"med.grpA.med.grpB"]) >= 1.5), "med.grpA.med.grpB"]),
			index.return=TRUE, decreasing=TRUE)
	
		gene.hits <- c(gene.hits.res,gene.hits.sens)
		
		gene.hits.sort <- c(gene.hits.res[gene.hits.res.sort$ix],gene.hits.sens[gene.hits.sens.sort$ix])
		
		if(length(gene.hits) < 2){
			next
		}
		
		# Add sort cell lines by tissue (within mut/non-mut) - sorted cell lines by tissue in kinome data
		# Add colour side bar to show tissues
		heatmap(
			t(zscores[grpAB, gene.hits.sort]),
#			breaks=c(-100,seq(from=-3, to=3, by=0.1),100),
			breaks=c(-100,seq(from=-3, to=0, by=0.05),100),
#			col=bluered(62),
			col=colorpanel(62,"darkblue","white"),
			ColSideColors=c(rep("lightgrey",length(grpA)),rep("darkblue",length(grpB))),
			Colv=NA,
			Rowv=NA,
			scale="none",
			cexRow=0.5,
			cexCol=0.8,
			margins=c(8,5),
			main=paste("mutated gene: ", gene, sep="")
		)
	}
	dev.off()
}


#
# make_heatmaps should be re-written to use image() instead of heatmap()
# then we can set the width appropriately to make each block square
#

make_heatmaps <- function(
	results,
	zscores,
	mutations,
	all_variants,
	outfile="zscore_heatmaps_by_mutated_genes.pdf"
	){	
	pdf(file=outfile)
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
#		gene <- strsplit(
#			colnames(mutations)[i],
#			"_"
#			)[[1]][1]
		gene <- colnames(mutations)[i]
		grpB <- which(all_variants[,gene] == 0)
		grpAB <- c(grpA,grpB)
		gene.hits.sens <- results[which(results[,"marker"] == gene & as.numeric(results[,"mptest.p"]) <= 0.05 & as.numeric(results[,"spearman.r"]) < -0.2 & as.numeric(results[,"nA"] >= 3)), 2]
		gene.hits.res <- results[which(results[,"marker"] == gene & as.numeric(results[,"mptest.p"]) <= 0.05 & as.numeric(results[,"spearman.r"]) > 0.2 & as.numeric(results[,"nA"] >= 3)), 2]
	
	
		gene.hits.sens.sort <- sort(
			as.numeric(results[which(results[,"marker"] == gene &
			as.numeric(results[,"mptest.p"]) <= 0.05 &
#			as.numeric(results[,"med.grpA.med.grpB"]) < 0), "med.grpA.med.grpB"]),
			as.numeric(results[,"spearman.r"]) < -0.2), "spearman.r"]),
			index.return=TRUE, decreasing=TRUE)
	
		gene.hits.res.sort <- sort(
			as.numeric(results[which(results[,"marker"] == gene &
			as.numeric(results[,"mptest.p"]) <= 0.05 &
#			as.numeric(results[,"med.grpA.med.grpB"]) > 0), "med.grpA.med.grpB"]),
			as.numeric(results[,"spearman.r"]) > 0.2), "spearman.r"]),
#			index.return=TRUE, decreasing=TRUE)
			index.return=TRUE, decreasing=TRUE)
	
		gene.hits <- c(gene.hits.res,gene.hits.sens)
		
		gene.hits.sort <- c(gene.hits.res[gene.hits.res.sort$ix],gene.hits.sens[gene.hits.sens.sort$ix])
		
		if(length(gene.hits) < 2){
			next
		}
		
		# Add sort cell lines by tissue (within mut/non-mut) - sorted cell lines by tissue in kinome data
		# Add colour side bar to show tissues
		heatmap(
			t(zscores[grpAB, gene.hits.sort]),
#			breaks=c(-100,seq(from=-3, to=3, by=0.1),100),
			breaks=c(-100,seq(from=-3, to=0, by=0.05),100),
#			col=bluered(62),
			col=colorpanel(62,"darkblue","white"),
			ColSideColors=c(rep("lightgrey",length(grpA)),rep("darkblue",length(grpB))),
			Colv=NA,
			Rowv=NA,
			scale="none",
			cexRow=0.8,
			cexCol=0.8,
			margins=c(8,5),
			main=paste("mutated gene: ", gene, sep="")
		)
	}
	dev.off()
}



make_volcano_plot <- function(
	results,
	outfile="volcano_mpt_meddiff.pdf"
	){
	
	# set mptest.p to 0.0001 if it is zero 
	# otherwise we get inf when log10 transforming
	results[which(results[,"mptest.p"] == 0),"mptest.p"] <- 0.0001
	
	pdf(file=outfile, height=10, width=15)
	plot(
		NULL,
		NULL,
		xlim=c(
			min(results[,"med.grpA.med.grpB"]),
			max(results[,"med.grpA.med.grpB"])
			),
		ylim=c(0,4.5),
		pch=19,
		cex=0.5,
		col="#999999",
		xlab="difference in group median z-score (mutant - non-mutant group)",
		ylab="-log10 p-value",
		cex.axis=1.5,
		cex.lab=1.5
		)
	redblackgrgreen.cols <- colorRampPalette(c("red", "black", "green"), interpolate="spline")
	my.redblackgreen.cols <- redblackgrgreen.cols(61)
	nmr.zscore.cuts <- seq(from=-3, to=3, by=0.1)
	nmr.zscore.cuts <- c(-100, nmr.zscore.cuts, 100)
	i <- NULL
	for(i in 2:length(nmr.zscore.cuts)){
		rows2plot <- which(results[,"med.grpB"] <= nmr.zscore.cuts[i] & results[,"med.grpB"] > nmr.zscore.cuts[i-1])
		
		
		points(
			results[rows2plot,"med.grpA.med.grpB"],
#			results[rows2plot,"med.grpA"]
			-log10(results[rows2plot,"mptest.p"]),
			pch=19,
			cex=1,
			col=my.redblackgreen.cols[i]
			) 
	}
	lines(c(-1,-1), c(0,100), lty=2, col="orange")
	lines(c(-2,-2), c(0,100), lty=2, col="red")
	lines(c(-100,100),c(-log10(0.05),-log10(0.05)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.01),-log10(0.01)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.001),-log10(0.001)), lty=2, col="blue")
	text(-14, 1.35, "p = 0.05", col="blue", cex=1)
	text(-14, 2.05, "p = 0.01", col="blue", cex=1)
	text(-14, 3.05, "p = 0.001", col="blue", cex=1)
	text(-1.9, 4.55, "∆z = -2", col="red", cex=1, pos=2)
	text(-0.9, 4.55, "∆z = -1", col="orange", cex=1, pos=2)
	
	
	rows_to_label <- as.numeric(rownames(
		results[which(
			results[,"med.grpA.med.grpB"] <= -2 & results[,"mptest.p"] <= 0.01
			)
			,]))
	rows_to_label_ordered <- rows_to_label[order(results[rows_to_label,"mptest.p"])]
	label_ys <- seq(from=4.4, to=-100, by=-0.08)
	last_label_ys_used <- 1
	i <- NULL

#	for(i in 1:nrow(results)){
	for(i in rows_to_label_ordered){
		if(results[i,"med.grpA.med.grpB"] <= -2 & results[i,"mptest.p"] <= 0.01){

			dot_x <- results[i,"med.grpA.med.grpB"]
			dot_y <- (-log10(results[i,"mptest.p"]))
			label_x <- -4.3
			label_y <- label_ys[last_label_ys_used]						
			text(
				label_x,
				label_y,
				paste(
					results[i,"marker"],
					":",
					results[i,"target"],
					sep=""
					),
				pos=2,
				cex=0.8
				)
			lines(
				c(label_x, dot_x),
				c(label_y, dot_y)
				)
			last_label_ys_used <- last_label_ys_used + 1
		}
	}

dev.off()
}


make_wilcox_volcano <- function(
	results,
	outfile="volcano_wilcox_meddiff.pdf"
	){
	
	pdf(file=outfile, height=10, width=15)
	plot(
		NULL,
		NULL,
#		xlim=c(
#			min(results[,"med.grpA.med.grpB"], na.rm=TRUE),
#			max(results[,"med.grpA.med.grpB"], na.rm=TRUE)
#			),
		xlim=c(-1.0,0.5),
		ylim=c(0,max(-log10(results[,"wilcox.p"]))),
		pch=19,
		cex=0.5,
		col="#999999",
		xlab="Spearman rho",
		ylab="-log10 p-value",
		cex.axis=1.5,
		cex.lab=1.5
		)
	redblackgrgreen.cols <- colorRampPalette(c("red", "black", "green"), interpolate="spline")
	my.redblackgreen.cols <- redblackgrgreen.cols(61)
	nmr.zscore.cuts <- seq(from=-3, to=3, by=0.1)
	nmr.zscore.cuts <- c(-100, nmr.zscore.cuts, 100)
	i <- NULL
	for(i in 2:length(nmr.zscore.cuts)){
		rows2plot <- which(results[,"med.grpB"] <= nmr.zscore.cuts[i] & results[,"med.grpB"] > nmr.zscore.cuts[i-1])
		
		if(length(rows2plot) < 1){
			next
		}
		
		points(
			results[rows2plot,"spearman.r"],
#			results[rows2plot,"med.grpA.med.grpB"],
#			results[rows2plot,"med.grpA"]
			-log10(results[rows2plot,"wilcox.p"]),
			pch=19,
			cex=1,
			col=my.redblackgreen.cols[i]
			) 
	}
	lines(c(-0.3,-0.3), c(0,100), lty=2, col="orange")
	lines(c(-0.5,-0.5), c(0,100), lty=2, col="red")
	lines(c(-100,100),c(-log10(0.05),-log10(0.05)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.01),-log10(0.01)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.001),-log10(0.001)), lty=2, col="blue")
	text(-14, 1.35, "p = 0.05", col="blue", cex=1)
	text(-14, 2.05, "p = 0.01", col="blue", cex=1)
	text(-14, 3.05, "p = 0.001", col="blue", cex=1)
	text(-0.4, 4.1, "r = -0.5", col="red", cex=1, pos=2)
	text(-0.2, 4.1, "r = -0.3", col="orange", cex=1, pos=2)
	
	
	rows_to_label <- as.numeric(rownames(
		results[which(
			results[,"spearman.r"] <= -0.3 & results[,"wilcox.p"] <= 0.005
			)
			,]))
	rows_to_label_ordered <- rows_to_label[order(results[rows_to_label,"wilcox.p"])]
	label_ys <- seq(from=4, to=-100, by=-0.06)
	last_label_ys_used <- 1
	i <- NULL

#	for(i in 1:nrow(results)){
	for(i in rows_to_label_ordered){
		if(results[i,"spearman.r"] <= -0.3 & results[i,"wilcox.p"] <= 0.005){

			dot_x <- results[i,"spearman.r"]
			dot_y <- (-log10(results[i,"wilcox.p"]))
			label_x <- -0.5
			label_y <- label_ys[last_label_ys_used]						
			text(
				label_x,
				label_y,
				paste(
					results[i,"marker"],
					":",
					results[i,"target"],
					sep=""
					),
				pos=2,
				cex=0.6
				)
			lines(
				c(label_x, dot_x),
				c(label_y, dot_y)
				)
			last_label_ys_used <- last_label_ys_used + 1
		}
	}

dev.off()
}

# modified version to cope with within-tissue results
# scales and lable positions are different
make_wilcox_volcano2 <- function(
	results,
	outfile="volcano_wilcox_meddiff.pdf"
	){
	
	pdf(file=outfile, height=10, width=15)
	plot(
		NULL,
		NULL,
#		xlim=c(
#			min(results[,"med.grpA.med.grpB"], na.rm=TRUE),
#			max(results[,"med.grpA.med.grpB"], na.rm=TRUE)
#			),
		xlim=c(-2,1),
		ylim=c(0,max(-log10(results[,"wilcox.p"]))),
		pch=19,
		cex=0.5,
		col="#999999",
		xlab="Spearman rho",
		ylab="-log10 p-value",
		cex.axis=1.5,
		cex.lab=1.5
		)
	redblackgrgreen.cols <- colorRampPalette(c("red", "black", "green"), interpolate="spline")
	my.redblackgreen.cols <- redblackgrgreen.cols(61)
	nmr.zscore.cuts <- seq(from=-3, to=3, by=0.1)
	nmr.zscore.cuts <- c(-100, nmr.zscore.cuts, 100)
	i <- NULL
	for(i in 2:length(nmr.zscore.cuts)){
		rows2plot <- which(results[,"med.grpB"] <= nmr.zscore.cuts[i] & results[,"med.grpB"] > nmr.zscore.cuts[i-1])
		
		if(length(rows2plot) < 1){
			next
		}
		
		points(
			results[rows2plot,"spearman.r"],
#			results[rows2plot,"med.grpA.med.grpB"],
#			results[rows2plot,"med.grpA"]
			-log10(results[rows2plot,"wilcox.p"]),
			pch=19,
			cex=1,
			col=my.redblackgreen.cols[i]
			) 
	}
	lines(c(-0.3,-0.3), c(0,100), lty=2, col="orange")
	lines(c(-0.5,-0.5), c(0,100), lty=2, col="red")
	lines(c(-100,100),c(-log10(0.05),-log10(0.05)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.01),-log10(0.01)), lty=2, col="blue")
	lines(c(-100,100),c(-log10(0.001),-log10(0.001)), lty=2, col="blue")
	text(-14, 1.35, "p = 0.05", col="blue", cex=1)
	text(-14, 2.05, "p = 0.01", col="blue", cex=1)
	text(-14, 3.05, "p = 0.001", col="blue", cex=1)
	text(-0.4, 4.1, "r = -0.5", col="red", cex=1, pos=2)
	text(-0.2, 4.1, "r = -0.3", col="orange", cex=1, pos=2)
	
	
	rows_to_label <- as.numeric(rownames(
		results[which(
			results[,"spearman.r"] <= -0.3 & results[,"wilcox.p"] <= 0.005
			)
			,]))
	rows_to_label_ordered <- rows_to_label[order(results[rows_to_label,"wilcox.p"])]
	label_ys <- seq(from=3.5, to=-100, by=-0.06)
	last_label_ys_used <- 1
	i <- NULL

#	for(i in 1:nrow(results)){
	for(i in rows_to_label_ordered){
		if(results[i,"spearman.r"] <= -0.3 & results[i,"wilcox.p"] <= 0.005){

			dot_x <- results[i,"spearman.r"]
			dot_y <- (-log10(results[i,"wilcox.p"]))
			label_x <- -0.9
			label_y <- label_ys[last_label_ys_used]						
			text(
				label_x,
				label_y,
				paste(
					results[i,"marker"],
					":",
					results[i,"target"],
					sep=""
					),
				pos=2,
				cex=0.6
				)
			lines(
				c(label_x, dot_x),
				c(label_y, dot_y)
				)
			last_label_ys_used <- last_label_ys_used + 1
		}
	}

dev.off()
}



make_spearman_volcano <- function(
		results,
		outfile="volcano_spearman.pdf"
	){
	pdf(file=outfile, width=13, height=9)
	plot(
		as.numeric(results[,"spearman.r"]),
		-log10(as.numeric(results[,"spearman.p"])),
		pch=19,
		cex=1,
		xlab="Spearman correlation coefficient",
		ylab="-log10(p-value)"
		)
	lines(
		c(-0.3,-0.3),
		c(0,10),
		lty=2,col="red"
		)
	lines(
		c(0.3,0.3),
		c(0,10),
		lty=2,col="green"
		)
	lines(
		c(-1,1),
		c(2,2),
		lty=2,col="gray"
		)
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results[i,17])){
			next;
		}
		if(as.numeric(results[i,15]) <= -0.3 & as.numeric(results[i,16]) <= 0.005){
				text(
					as.numeric(results[i,"spearman.r"]),
					-log10(as.numeric(results[i,"spearman.p"])),
					paste(results[i,1],results[i,2], sep=":"),
					cex=1,
					pos=4
					)
			}
	}
	dev.off()
	1
}


# function from http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# to draw a scale bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=round(seq(min, max, len=nticks), 2), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, cex.main=2)
    axis(2, ticks, las=1, cex=2, cex.axis=2)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}




# This has been superceeded by plot_marker_target4
# It also won't work until options are added to feed in the data for the mutations and z-scores...
plot_marker_target3 <- function(marker, target){
	grpA <- which(wtsi.genes.func.cmn[,marker] > 0)
	grpB <- which(wtsi.excl.cmn[,marker] == 0)
	
	target.zscores <- 
		rbind(
			cbind(
				paste(marker, "mutant"),
				na.omit(kinome.cmn[grpA, target])
			),
			cbind(
				paste(marker, "non-mutant"),
				na.omit(kinome.cmn[grpB, target])
			)
		)

	png(file=paste(marker,"_mutants_by_", target,"_zscores_131218.png", sep=""),
		width=600,
		height=480,
		pointsize=18
		)
	stripchart(
		as.numeric(target.zscores[,2])~target.zscores[,1],
		pch=19,
		col="grey",
		vertical=TRUE,
		method="jitter",
		jitter=0.01,
		ylab=paste(target, "z-score"),
		xlim=c(0.5,2.5),
		cex.axis=1.2,
		cex.lab=1.2
		)
	i <- NULL
	labels_and_ys <- data.frame(c(NA,NA),stringsAsFactors=FALSE) # store the mutations details and y-positions for each label
	labels_and_ys_row <- 0
	mutants <- NULL
	for(i in seq(1:ncol(wtsi.cmn))){
		gene <- strsplit(
			colnames(wtsi.cmn)[i],
			"_"
			)[[1]][1]
		if(gene == marker){
			mutants <- which(wtsi.cmn[,i] > 0)
			
			j <- NULL
			for(j in 1:length(mutants)){
				mutation <- strsplit(colnames(wtsi.cmn)[i],"_")[[1]][2]
				if(mutation == "p.?"){
					next
				}
				if(mutation == "NA"){
					next
				}
				labels_and_ys_row <- labels_and_ys_row + 1
				labels_and_ys[labels_and_ys_row,1] <- mutation
				labels_and_ys[labels_and_ys_row,2] <- as.numeric(kinome.cmn[mutants[j],target])
			}
		}
	}
	# get rid of any cell lines with NAs for z-scores...
	labels_and_ys <- na.omit(labels_and_ys)
	
	
	last_y <- 1000 # set the last y-pos to a very big number and update in loop
	# the value of 23.8 below was found empirically by printing out a plot and 
	# measuring. A capital letter is 0.5 cm and the plot height is 11.9 cm.
	# Therefore there are 23.8 capitals that could fit into the height of the plot
	# This was increased to 24 to fit a bit better
	min_lable_dist <- as.numeric((max(na.omit(kinome.cmn[,target])) - min(na.omit(kinome.cmn[,target]))) / 25)


	# need to sort mutants by z-scores...
	labels_and_ys.sort.order <- sort(labels_and_ys[,2], decreasing=TRUE, index.return=TRUE)$ix
#	mutants <- mutants[mutant.sort.order] # order mutants by decreasing z-value

	k <- NULL
	for(k in labels_and_ys.sort.order){
		text_y <- NULL # this is where the text will be drawn
		if(labels_and_ys[k,2] > (as.numeric(last_y) - as.numeric(min_lable_dist)) ){ 
			text_y <- labels_and_ys[k,2] - (min_lable_dist - (as.numeric(last_y) - as.numeric(labels_and_ys[k,2])))
		}else{
			text_y <- as.numeric(labels_and_ys[k,2])
		}
		text(
			1.05,
			text_y,
			labels_and_ys[k,1],
			pos=4,
			cex=0.7
			)
	
		lines(
			c(1.02,1.1),
			c(labels_and_ys[k,2], text_y),
			lwd=2
			)
		last_y <- text_y		
	}

	dev.off()	
}




# The function below is broken unless we also import the complete mutation detais...
#
#
# add more info to the plots (line for sensitivity threshold, p-values etc)
# also, check if target z-scores exist - if NA, don't plot the mutation information - No!
plot_marker_target4 <- function(
	marker,
	target,
	threshold,
	pval,
	medA,
	medB,
	zscores,
	mutations,
	all_variants,
	outfile="zscore_stripcharts_by_mutated_genes.pdf"
	){

#	just make the image and rely on the caller to open the PDF	
#	pdf(file=outfile")
	
	grpA <- which(mutations[,marker] > 0)
	grpB <- which(all_variants[,marker] == 0)
	
	target.zscores <- 
		rbind(
			cbind(
				paste(marker, "mutant"),
				na.omit(zscores[grpA, target])
			),
			cbind(
				paste(marker, "non-mutant"),
				na.omit(zscores[grpB, target])
			)
		)
	
	pval.text <- NULL
	if(round(pval, 3) == 0){
		pval.text <- "<0.000"
	}else{
		pval.text <- paste("=", round(pval, 3), sep="")
	}
	
#	png(file=paste(marker,"_mutants_by_", target,"_zscores_131218.png", sep=""),
#		width=600,
#		height=480,
#		pointsize=18
#		)
	stripchart(
		as.numeric(target.zscores[,2])~target.zscores[,1],
		pch=19,
		col="grey",
		vertical=TRUE,
		method="jitter",
		jitter=0.01,
		ylab=paste(target, "z-score"),
		xlim=c(0.5,2.5),
		cex.axis=1.2,
		cex.lab=1.2,
		main=paste(
			"Mutation: ",
			marker,
			", Group medians: ",
			round(medA, 2),
			" (mutant), ",
			round(medB, 2),
			" (non-mutant)\n",
			"(Wilcoxon Rank Sum test: p",
			pval.text,
			")",
			sep=""
			)
		)
	lines(
		c(0.5,2.5),
		c(threshold, threshold),
		col="red",
		lty=2
		)
	lines(
		c(0.875, 1.125),
		c(medA, medA),
		col="blue",
		lty=2
		)
	lines(
		c(1.875, 2.125),
		c(medB, medB),
		col="blue",
		lty=2
		)
	
	
	
	i <- NULL
	labels_and_ys <- data.frame(c(NA,NA),stringsAsFactors=FALSE) # store the mutations details and y-positions for each label
	labels_and_ys_row <- 0
	mutants <- NULL
	for(i in seq(1:ncol(wtsi.cmn))){
		gene <- strsplit(
			colnames(wtsi.cmn)[i],
			"_"
			)[[1]][1]
		if(gene == marker){
			mutants <- which(wtsi.cmn[,i] > 0)
			
			j <- NULL
			for(j in 1:length(mutants)){
				mutation <- strsplit(colnames(wtsi.cmn)[i],"_")[[1]][2]
				if(mutation == "p.?"){
					next
				}
				if(mutation == "NA"){
					next
				}
				labels_and_ys_row <- labels_and_ys_row + 1
				labels_and_ys[labels_and_ys_row,1] <- mutation
				labels_and_ys[labels_and_ys_row,2] <- as.numeric(zscores[mutants[j],target])
			}
		}
	}
	# get rid of any cell lines with NAs for z-scores...
	labels_and_ys <- na.omit(labels_and_ys)
	
	
	last_y <- 1000 # set the last y-pos to a very big number and update in loop
	# the value of 23.8 below was found empirically by printing out a plot and 
	# measuring. A capital letter is 0.5 cm and the plot height is 11.9 cm.
	# Therefore there are 23.8 capitals that could fit into the height of the plot
	# This was increased to 24 to fit a bit better
	min_lable_dist <- as.numeric((max(na.omit(zscores[,target])) - min(na.omit(zscores[,target]))) / 25)


	# need to sort mutants by z-scores...
	labels_and_ys.sort.order <- sort(labels_and_ys[,2], decreasing=TRUE, index.return=TRUE)$ix
#	mutants <- mutants[mutant.sort.order] # order mutants by decreasing z-value

	k <- NULL
	for(k in labels_and_ys.sort.order){
		text_y <- NULL # this is where the text will be drawn
		if(labels_and_ys[k,2] > (as.numeric(last_y) - as.numeric(min_lable_dist)) ){ 
			text_y <- labels_and_ys[k,2] - (min_lable_dist - (as.numeric(last_y) - as.numeric(labels_and_ys[k,2])))
		}else{
			text_y <- as.numeric(labels_and_ys[k,2])
		}
		text(
			1.05,
			text_y,
			labels_and_ys[k,1],
			pos=4,
			cex=0.7
			)
	
		lines(
			c(1.02,1.1),
			c(labels_and_ys[k,2], text_y),
			lwd=2
			)
		last_y <- text_y		
	}

#	We are expecting the caller to open the graphics device
#	dev.off()	
}





run_univariate_tests <- function(
	zscores,
	mutations,
	all_variants,
	sensitivity_thresholds=NULL,
	nperms=1000000,
	alt="less"
	){
	
	
	zscores <- as.matrix(zscores)
	mutations <- as.matrix(mutations)
	all_variants <- as.matrix(all_variants)
	
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		
		#gene <- strsplit(
		#	colnames(mutations)[i],
		#	"_"
		#	)[[1]][1]
		gene <- colnames(mutations)[i]
		
		
		# grpB includes cell lines with no reported mutations at all
		# in gene...
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3 as we are never going to
		# consider anything based on n=2
		if(length(grpA) < 3 | length(grpB) < 3){
			next
		}
		

# this code changed because it does not exclude
# uncertain observations
#		# this is used for spearman correlation...
#		mut.status <- rep(0,nrow(zscores))
#		mut.status[which(mutations[,i] > 0)] <- 1

		# this is used for spearman correlation...
		mut.status <- rep(NA,nrow(zscores))
		mut.status[which(mutations[,i] == 1)] <- 1
		mut.status[which(all_variants [,i] == 0)] <- 0

		
		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){
			
			# skip if we have no viability measurements
			# in one or other group
			if(length(na.omit(zscores[grpA,j])) < 3){
				next
			}
			if(length(na.omit(zscores[grpB,j])) < 3){
				next
			}
			
			# calc the median permutation test p-value
			real.med.diff <- median(na.omit(zscores[grpA,j])) - median(na.omit(zscores[grpB,j]))
			
			# permute the groups nperms times, sampling with group sizes equal to grpA
			sample.size <- length(grpA)

			mpt.pval <- "NA"

# begin uncommenting code for MPTest

			k <- NULL
			permuted.med.diffs <- NULL
			grpAB <- c(grpA,grpB) # join up the actual cell lines we used so we can sample them.

			for(k in 1:1000){
				index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
				permuted.grpA <- grpAB[index]
				permuted.grpB <- grpAB[-index]
				permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
			}
			mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 1000		

			if(mpt.pval < 0.010){
				for(k in 1001:10000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 10000
			}

			if(mpt.pval < 0.0010){
				for(k in 10001:100000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 100000
			}

			if(mpt.pval < 0.00010){
				for(k in 100001:1000000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 1000000
			}

# end uncommenting code for MPTest
			
			wilcox.p <- NA
			try(
				test <- wilcox.test(
					zscores[grpA,j],
					zscores[grpB,j],
					alternative=alt
				)
			)
			wilcox.p <- test$p.value
			
			# get the Spearman r value as an effect size estimate
			spearman <- NULL
			try(
				spearman <- cor.test(zscores[,j], mut.status, method="spearman", use="complete.obs", alternative=alt)
			)
			
			this_threshold <- 0
			if(!is.null(sensitivity_thresholds)){
				this_threshold <- sensitivity_thresholds[colnames(zscores)[j]]
			}
			
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			mad.grpA <- mad(zscores[grpA,j], na.rm=TRUE)
			mad.grpB <- mad(zscores[grpB,j], na.rm=TRUE)
			countA.sens <- length(which(zscores[grpA,j] <= this_threshold))
			countB.sens <- length(which(zscores[grpB,j] <= this_threshold))
#			countA.sens <- length(which(zscores[grpA,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
#			countB.sens <- length(which(zscores[grpB,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			nMin <- min(nA, nB)
			pcnt.grpA.sens <- 100 * countA.sens / nA
			pcnt.grpB.sens <- 100 * countB.sens / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			spearman.r <- spearman$estimate
			spearman.p <- spearman$p.value
			
			# Output the result if min sample size is 2 or more
			if(nMin > 1){
				results <- rbind(
					results,
					c(
						marker,
						target,
						nA,
						nB,
						this_threshold,
						med.grpA,
						med.grpB,
						med.diff,
						countA.sens,
						countB.sens,
						pcnt.grpA.sens,
						pcnt.grpB.sens,
						min.grpA,
						min.grpB,
						spearman.r,
						spearman.p,
						wilcox.p,
						mpt.pval
					)
				)
			}
		}
	}
	
	if(is.null(nrow(results))){
		return(NULL)
	}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
		"percent.grpA.sens",
		"percent.grpB.sens",
		"min.grpA",
		"min.grpB",
		"spearman.r",
		"spearman.p",
		"wilcox.p",
		"mptest.p"
	)
	
	return(results)
	
} # end run_univariate_tests


#
# Test kinome targets with CNA - modification of run_univariate_tests
# only test self with self
#



test_kinome_with_self_cnv <- function(
	zscores,
	mutations,
	all_variants,
	sensitivity_thresholds,
	nperms=1000
	){
	
	
	zscores <- as.matrix(zscores)
	mutations <- as.matrix(mutations)
	all_variants <- as.matrix(all_variants)
	
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		
		#gene <- strsplit(
		#	colnames(mutations)[i],
		#	"_"
		#	)[[1]][1]
		gene <- colnames(mutations)[i]
		
		
		# grpB includes cell lines with no reported mutations at all
		# in gene...
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3 as we are never going to
		# consider anything based on n=2
		if(length(grpA) < 3 | length(grpB) < 3){
			next
		}




# this code changed because it does not exclude
# uncertain observations
#		# this is used for spearman correlation...
#		mut.status <- rep(0,nrow(zscores))
#		mut.status[which(mutations[,i] > 0)] <- 1

		# this is used for spearman correlation...
		mut.status <- rep(NA,nrow(zscores))
		mut.status[which(mutations[,i] == 1)] <- 1
		mut.status[which(all_variants [,i] == 0)] <- 0

		j <- which(colnames(zscores) == gene)
			
			# calc the median permutation test p-value
			abs.med.diff <- abs(median(na.omit(zscores[grpA,j])) - median(na.omit(zscores[grpB,j])))
			
			# permute the groups nperms times, sampling with group sizes equal to grpA
			sample.size <- length(grpA)
#			k <- NULL
#			permuted.abs.med.diffs <- NULL
#			grpAB <- c(grpA,grpB) # join up the actual cell lines we used so we can sample them.
#			
#			for(k in 1:nperms){
#				index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
#				permuted.grpA <- grpAB[index]
#				permuted.grpB <- grpAB[-index]
#				permuted.abs.med.diffs[k] <- abs(median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j])))
#			}
#			mpt.pval <- length(which(permuted.abs.med.diffs >= abs.med.diff)) / nperms
			
			mpt.pval <- "NA"
			
			wilcox.p <- NA
			try(
				test <- wilcox.test(
					zscores[grpA,j],
					zscores[grpB,j],
					alternative="less"
				)
			)
			wilcox.p <- test$p.value
			
			# get the Spearman r value as an effect size estimate
			spearman <- NULL
			spearman <- cor.test(zscores[,j], mut.status, method="spearman", use="complete.obs")
			
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			mad.grpA <- mad(zscores[grpA,j], na.rm=TRUE)
			mad.grpB <- mad(zscores[grpB,j], na.rm=TRUE)
			countA.sens <- length(which(zscores[grpA,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
			countB.sens <- length(which(zscores[grpB,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			nMin <- min(nA, nB)
			pcnt.grpA.sens <- 100 * countA.sens / nA
			pcnt.grpB.sens <- 100 * countB.sens / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			spearman.r <- spearman$estimate
			spearman.p <- spearman$p.value
			
			# Output the result if min sample size is 2 or more
			if(nMin > 1){
				results <- rbind(
					results,
					c(
						marker,
						target,
						nA,
						nB,
						sensitivity_thresholds[colnames(zscores)[j]],
						med.grpA,
						med.grpB,
						med.diff,
						countA.sens,
						countB.sens,
						pcnt.grpA.sens,
						pcnt.grpB.sens,
						min.grpA,
						min.grpB,
						spearman.r,
						spearman.p,
						wilcox.p,
						mpt.pval
					)
				)
			}
	}
	
	if(is.null(nrow(results))){
		return(NULL)
	}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
		"percent.grpA.sens",
		"percent.grpB.sens",
		"min.grpA",
		"min.grpB",
		"spearman.r",
		"spearman.p",
		"wilcox.p",
		"mptest.p"
	)
	
	return(results)
	
} 


###############

#
# version of run_univariate_tests() used for validation screens (t-test)
#

run_validation_tests <- function(
	zscores,
	mutations,
	all_variants,
	sensitivity_thresholds,
	nperms=1000000
	){
	
	
	zscores <- as.matrix(zscores)
	mutations <- as.matrix(mutations)
	all_variants <- as.matrix(all_variants)
	
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		grpA <- which(mutations[,i] > 0)
		
		#gene <- strsplit(
		#	colnames(mutations)[i],
		#	"_"
		#	)[[1]][1]
		gene <- colnames(mutations)[i]
		
		
		# grpB includes cell lines with no reported mutations at all
		# in gene...
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3 as we are never going to
		# consider anything based on n=2
		if(length(grpA) < 3 | length(grpB) < 3){
			next
		}

		# this is used for spearman correlation...
		mut.status <- rep(NA,nrow(zscores))
		mut.status[which(mutations[,i] == 1)] <- 1
		mut.status[which(all_variants [,i] == 0)] <- 0

		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){
			
			# calc the median permutation test p-value
			real.med.diff <- median(na.omit(zscores[grpA,j])) - median(na.omit(zscores[grpB,j]))
			
			# permute the groups nperms times, sampling with group sizes equal to grpA
			sample.size <- length(grpA)
			k <- NULL
			permuted.med.diffs <- NULL
			grpAB <- c(grpA,grpB) # join up the actual cell lines we used so we can sample them.
			
			# modified to run a one-sided test. Only asking if the permuted groups median is
			# ≤ the actual median
			for(k in 1:nperms){
				index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
				permuted.grpA <- grpAB[index]
				permuted.grpB <- grpAB[-index]
				permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
			}
			mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / nperms
			
			wilcox.p <- NA
			try(
				test <- t.test(
					zscores[grpA,j],
					zscores[grpB,j],
					alternative="less"
				)
			)
			wilcox.p <- test$p.value
			
			# get the Spearman r value as an effect size estimate
			spearman <- NULL
			try(
				spearman <- cor.test(zscores[,j], mut.status, method="spearman", use="complete.obs", alternative="less")
			)
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			mad.grpA <- mad(zscores[grpA,j], na.rm=TRUE)
			mad.grpB <- mad(zscores[grpB,j], na.rm=TRUE)
			countA.sens <- length(which(zscores[grpA,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
			countB.sens <- length(which(zscores[grpB,j] <= sensitivity_thresholds[colnames(zscores)[j]]))
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			nMin <- min(nA, nB)
			pcnt.grpA.sens <- 100 * countA.sens / nA
			pcnt.grpB.sens <- 100 * countB.sens / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			spearman.r <- spearman$estimate
			spearman.p <- spearman$p.value
			
			# Output the result if min sample size is 2 or more
			if(nMin > 1){
				results <- rbind(
					results,
					c(
						marker,
						target,
						nA,
						nB,
						sensitivity_thresholds[colnames(zscores)[j]],
						med.grpA,
						med.grpB,
						med.diff,
						countA.sens,
						countB.sens,
						pcnt.grpA.sens,
						pcnt.grpB.sens,
						min.grpA,
						min.grpB,
						spearman.r,
						spearman.p,
						wilcox.p,
						mpt.pval
					)
				)
			}
		}
	}
	
	if(is.null(nrow(results))){
		return(NULL)
	}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
		"percent.grpA.sens",
		"percent.grpB.sens",
		"min.grpA",
		"min.grpB",
		"spearman.r",
		"spearman.p",
		"ttest.p",
		"mptest.p"
	)
	
	return(results)
	
} # end run_univariate_tests

count_z_lt_neg2 <- function(x){
		length(which(x < -2))
}

#
# Run the full set of analyses used for Intercell
#

run_intercell <- function(
	x,
	qnorm=FALSE,
	analysis_id="v1_yymmdd"
	){
	# Run a standard analysis for a given list of data tables
	
	# load modules used in parts of the analysis if not already done
	require("preprocessCore")
	require(gplots)
	
	# first decide if we are going to be working with Z-scores or 
	# quantile normalised Z-scores
	zscores <- x$rnai
	sensitivity_thresholds <- x$rnai_iqr_thresholds
	if(qnorm == TRUE){
		zscores <- x$rnai_qn
		sensitivity_thresholds <- x$rnai_qn_iqr_thresholds
	}

	
	# define strong dependencies irrespective of genotype or histotype
	target_z_lt_neg2_counts <- apply(zscores, 2, count_z_lt_neg2)
	pdf(
		file=paste(
			"freq_plot_targets_z_lessthat_neg2_",
			analysis_id,
			".pdf",
			sep=""
			),
		width=4.5,
		height=4
		)
	plot(
		sort(target_z_lt_neg2_counts, decreasing=TRUE),
		pch=19,
		col=rgb(0,0,0,0.5),
		xlab="siRNA targets",
		ylab="number of dependent cell lines",
		las=1,
		log="x"
		)
	# line at 10% (3.8)
	x10pc_line_y <- nrow(zscores) * 0.1
	x30pc_line_y <- nrow(zscores) * 0.3
	x50pc_line_y <- nrow(zscores) * 0.5
	x70pc_line_y <- nrow(zscores) * 0.7
	x90pc_line_y <- nrow(zscores) * 0.9
	abline(
		x10pc_line_y,
		0,
		lty=2
		)
	text(
		ncol(zscores) - 26,
		x10pc_line_y - 1,
		"10%",
		pos=3
		)
	
	# line at 30% (11.4)
	abline(
		x30pc_line_y,
		0,
		lty=2
		)
	text(
		ncol(zscores) - 26,
		x30pc_line_y - 1,
		"30%",
		pos=3
		)
	
	# line at 50% (19)
	abline(
		x50pc_line_y,
		0,
		lty=2
		)
	text(
		ncol(zscores) - 26,
		x50pc_line_y - 1,
		"50%",
		pos=3
		)
	
	# line at 70% (26.6)
	abline(
		x70pc_line_y,
		0,
		lty=2
		)
	text(
		ncol(zscores) - 26,
		x70pc_line_y - 1,
		"70%",
		pos=3
		)
	# line at 90% (26.6)
	abline(
		x90pc_line_y,
		0,
		lty=2
		)
	text(
		ncol(zscores) - 26,
		x90pc_line_y - 1,
		"90%",
		pos=3
		)
		
	# dummy plot
	plot(NULL,NULL,xlim=c(1,0), ylim=c(0,1))
	dev.off()
	
	num_cell_lines_sensitive_counts <- NULL
	
	num_cell_lines_sensitive_counts <- rbind(
		c(
			"1 or more",
			length(which(target_z_lt_neg2_counts > 0))
			),
		c(
			"5 or more",
			length(which(target_z_lt_neg2_counts > 4))
			),
		c(
			"10 or more",
			length(which(target_z_lt_neg2_counts > 9))
			),
		c(
			"50% or more",
			length(which(target_z_lt_neg2_counts > (nrow(zscores)*0.5)))
			),
		c(
			"90% or more",
			length(which(target_z_lt_neg2_counts > (nrow(zscores)*0.9)))
			)
		)

	colnames(num_cell_lines_sensitive_counts) <- c(
		"group",
		"number_of_cell_lines"
		)
	write.table(
		num_cell_lines_sensitive_counts,
		file=paste(
			"freq_table_targets_z_lessthat_neg2_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE,
		col.names=TRUE
		)
	
	# 
	# look for dependencies associated with histotype
	#
	
	
	# set up colour scheme for tissue types
	legend_pretty_tissues = c(
		"Osteosarcoma",
		"Breast",
		"Lung",
		"Head & Neck",
		"Pancreas",
		"Cervical",
		"Ovary",
		"Esophagus",
		"Endometrium",
		"CNS"
		)
	legend_actual_tissues = c(
		"BONE",
		"BREAST",
		"LUNG",
		"HEADNECK",
		"PANCREAS",
		"CERVICAL",
		"OVARY",
		"OESOPHAGUS",
		"ENDOMETRIUM",
		"CENTRAL_NERVOUS_SYSTEM"
		)
	legend_col=c(
		"yellow",
		"deeppink",
		"darkgrey",
		"firebrick4",
		"purple",
		"blue",
		"cadetblue",
		"green",
		"orange",
		"darkgoldenrod4"
		)
	names(legend_col) <- legend_actual_tissues
	
	tissues_used <- which(legend_actual_tissues %in% colnames(x$tissues))
	
	legend_pretty_tissues <- legend_pretty_tissues[tissues_used]
	legend_actual_tissues <- legend_actual_tissues[tissues_used]
	legend_col <- legend_col[tissues_used]
	
	x_tissues <- x
	
	x_tissues$func_muts <- x$tissues	
	x_tissues$all_muts <- x$tissues	
	x_tissues$mut_classes <- x$tissues	
	
	uv_results_x_tissues <- run_univariate_tests(
		zscores=zscores,
		mutations=x_tissues$func_muts,
		all_variants=x_tissues$all_muts,
		sensitivity_thresholds=sensitivity_thresholds
		)
	
	write.table(
		uv_results_x_tissues,
		file=paste(
			"dependencies_associated_with_histotype_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE
	)

	uv_results_x_tissues <- read.table(
		file=paste(
			"dependencies_associated_with_histotype_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		header=TRUE,
		stringsAsFactors=FALSE
	)


	#
	# dependencies associated with driver gene mutation status
	#
	
	# write the mutation tables for reference
	write.table(
		x$func_muts,
		file=paste(
			"functional_mutation_dataset_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE
	)
	write.table(
		x$all_muts,
		file=paste(
			"uncertain_function_mutation_dataset_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE
	)
	
	uv_results_x_combuts_in_all_histotypes <- run_univariate_tests(
		zscores=zscores,
		mutations=x$func_muts,
		all_variants=x$all_muts,
		sensitivity_thresholds=sensitivity_thresholds
		)

	write.table(
		uv_results_x_combuts_in_all_histotypes,
		file=paste(
			"dependencies_associated_with_driver_mutations_in_all_histotypes_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE
		)

	uv_results_x_combuts_in_all_histotypes <- read.table(
		file=paste(
			"dependencies_associated_with_driver_mutations_in_all_histotypes_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		header=TRUE,
		stringsAsFactors=FALSE
		)

	#
	# dependencies associated with driver gene mutations in each histotype
	#
	
	tissue_types <- colnames(x$tissues)
	
	uv_results_bytissue <- NULL
	tissue <- NULL
	for(tissue in tissue_types){
		
		cellline_count <- sum(
			x$tissues[,tissue]
			)
		
		if(cellline_count < 5){
			next
		}
		
		tissue_rows <- which(
			x$tissues[,tissue] == 1
			)
		
		temp_results <- NULL
		temp_results <- run_univariate_tests(
			zscores=zscores[tissue_rows,],
			mutations=x$func_muts[tissue_rows,],
			all_variants=x$all_muts[tissue_rows,],
			sensitivity_thresholds=sensitivity_thresholds
			)
		
		if(is.null(nrow(temp_results))){
			print(paste("Skipping ", tissue, " - no results", sep=""))
			next
		}
		
		temp_results <- cbind(
			temp_results,
			rep(tissue, times=nrow(temp_results))
			)
		
		uv_results_bytissue <- rbind(
			uv_results_bytissue,
			temp_results
			)
	}
	colnames(
		uv_results_bytissue
		)[ncol(
			uv_results_bytissue
			)
			] <- "tissue"
	
	
	write.table(
		uv_results_bytissue,
		file=paste(
			"dependencies_associated_with_driver_mutations_in_individual_histotypes_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		quote=FALSE,
		row.names=FALSE
		)

	uv_results_bytissue <- read.table(
		file=paste(
			"dependencies_associated_with_driver_mutations_in_individual_histotypes_",
			analysis_id,
			".txt",
			sep=""
			),
		sep="\t",
		header=TRUE,
		stringsAsFactors=FALSE
		)
	
	#
	# make box plots of significant and all associations
	#
	
	
	# Make mini skinny boxplots for combined histotype assocs
	make_mini_box_dot_plots4_col_by_tissue(
		results=as.data.frame(
			uv_results_x_tissues
			),
		zscores=zscores,
		mutation.classes=x_tissues$mut_classes,
		mutations=x_tissues$func_muts,
		exclusions=x_tissues$all_muts,
		tissues=x_tissues$tissues,
		filename=paste(
			"boxplots_dependencies_associated_with_histotype_",
			analysis_id,
			".pdf",
			sep=""
			),
		tissue_pretty_names=legend_pretty_tissues,
		tissue_actual_names=legend_actual_tissues,
		tissue_cols=legend_col
		)


	make_mini_box_dot_plots4_col_by_tissue(
		results=as.data.frame(
			uv_results_x_combuts_in_all_histotypes
			),
		zscores=zscores,
		mutation.classes=x$mut_classes,
		mutations=x$func_muts,
		exclusions=x$all_muts,
		tissues=x$tissues,
		filename=paste(
			"boxplots_dependencies_associated_with_driver_mutations_in_all_histotypes_",
			analysis_id,
			".pdf",
			sep=""
			),
		tissue_pretty_names=legend_pretty_tissues,
		tissue_actual_names=legend_actual_tissues,
		tissue_cols=legend_col
		)

	
	# Make mini skinny boxplots for histotype-specific assocs
	tissues <- levels(as.factor(uv_results_bytissue$tissue))
	for(tissue in tissues){
		results_tissue_rows <- which(uv_results_bytissue$tissue == tissue)
		data_tissue_rows <- which(x$tissues[,tissue] == 1)
		make_mini_box_dot_plots4_col_by_tissue(
			results=as.data.frame(
				uv_results_bytissue[results_tissue_rows,]
				),
			zscores=zscores[data_tissue_rows,],
			mutation.classes=x$mut_classes[data_tissue_rows,],
			mutations=x$func_muts[data_tissue_rows,],
			exclusions=x$all_muts[data_tissue_rows,],
			tissues=x$tissues[data_tissue_rows,],
			filename=paste(
				"boxplots_dependencies_associated_with_driver_mutations_in_",
				tissue,
				"_",
				analysis_id,
				".pdf",
				sep=""
				),
			tissue_pretty_names=legend_pretty_tissues,
			tissue_actual_names=legend_actual_tissues,
			tissue_cols=legend_col
		)
	}
	
	
	# cheese and pineapples
	
	# driver mutations barplot

}




















