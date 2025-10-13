library(nebula)
library(Matrix)
library(Seurat)
library(tidyverse)

start_time = proc.time()

nebula_output_dir <- "nebula_output_dir/"
dir.create(nebula_output_dir)

unique_cell_types <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial",
		"Epicardial", "Fibroblast", "LEC", "Lymphoid", "Mast", "Myeloid",
		"Neuronal", "Pericyte", "vSMC")

for (cell_type_val in unique_cell_types) {

	print(cell_type_val)
	flush.console()
	
	seurat_obj_path <- paste0("Seurat_obj_dir/Seurat_", cell_type_val, "/seurat_object.rds")
	subset_obj <- readRDS(seurat_obj_path)
	
	# get log10 UMI
	#subset_obj@metadata$log_UMI <- log10(seurat_obj@meta.data$nCount_RNA)

	# convert Seurat to Nebula
	seuratdata <- scToNeb(obj = subset_obj, 
                      assay = "RNA", 
                      id = "donor_id", 
                      pred = c("sex", "age_group", "tech_plus_study"), 
                      offset="nCount_RNA")
	
	#print(seuratdata)

    	# make sure young is the reference for age group
	# female is the reference for sex
    	seuratdata$pred$age_group <- relevel(as.factor(seuratdata$pred$age_group), ref = "young")
    	seuratdata$pred$sex <- relevel(as.factor(seuratdata$pred$sex), ref = "female")
	#seuratdata$pred$disease_binary <- relevel(as.factor(seuratdata$pred$disease_binary), ref = "N")
    
    	#df = model.matrix(~ sex + age_group + tech_plus_study + disease_binary, data = seuratdata$pred)
	df = model.matrix(~ sex + age_group + tech_plus_study, data = seuratdata$pred)

	# order the cells by donor
	data_g = group_cell(count=seuratdata$count,
                    id=seuratdata$id,
                    pred=df,
                    offset=seuratdata$offset)

	print("Running nebula model...")
	flush.console()

	# run nebula
	re = nebula(data_g$count,
            data_g$id,
            pred=df,
            offset=data_g$offset,
            ncore=20)

	# free RAM periodically
	gc()

	# save the RDS 
	saveRDS(object = re, 
		file = paste0(nebula_output_dir, cell_type_val, ".rds"))
}

end_time = proc.time()
elapsed_time = end_time - start_time
print(paste0("Script complete! Elapsed time in seconds is ", elapsed_time))
