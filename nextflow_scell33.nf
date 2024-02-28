// Define input parameters
params {
    work_dir = "/Singlecell_Project33"
    script_dir = "$work_dir/nextflow_scell33"
    input_data_dir = "$script_dir/input_data"
    input_files = ["GSM5764336_filtered_feature_bc_matrix_ABU8.h5", "GSM5764394_filtered_feature_bc_matrix_CS161.h5", ...]
}

// Process to read 10x data
process read_10x {
    input:
    file input_file from input_files

    output:
    file "${input_file.baseName}_seurat_obj.rds" into seurat_obj_file

    script:
    """
    library(Seurat)
    library(Matrix)
    
    GSM_data <- Seurat::Read10X_h5("$input_data_dir/$input_file", use.names = TRUE, unique.features = TRUE)
    saveRDS(GSM_data, "${input_file.baseName}_seurat_obj.rds")
    """
}

// Process to preprocess data
process preprocess_data {
    input:
    file seurat_obj_file from read_10x

    output:
    file "PCA_plot.rds" into pca_plot_file,
    file "PP_PN_cluster_annotate.rds" into pp_pn_cluster_annotate_file,
    file "UMAP_plot.rds" into umap_plot_file,
    file "PTFgeneexpression_dotplots.rds" into ptf_geneexpression_dotplots_file,
    file "PRfeat.rds" into pr_feat_file,
    file "PNfeat.rds" into pn_feat_file,
    file "seurat_obj.rds" into final_seurat_obj_file

    script:
    """
    library(Seurat)
    
    # Load Seurat object
    seurat_obj <- readRDS("$seurat_obj_file")
    
    # Preprocessing data
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 3 & nFeature_RNA < 350)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    
    # Perform PCA
    seurat_obj <- RunPCA(seurat_obj, npcs = 20)
    
    # Find neighbors and clusters
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    
    # Gene lists
    your_gene_list_PR <- c("PTPRU", "PTPRF", "PTPRC", "PTPRN", "PTPRG", "PTPRK", "PTPRZ1", "PTPRN2", "PTPRD", "PTPRE", "PTPRJ", "PTPRCAP", "PTPRO", "PTPRB", "PTPRR", "PTPRQ", "PTPRM", "PTPRA", "PTPRT", "PTPRS", "PTPRH")
    your_gene_list_PN <- c("PTPN22", "PTPN7", "PTPN14", "PTPN4", "PTPN18", "PTPN23", "PTPN13", "PTPN12", "PTPN3", "PTPN20A", "PTPN20B", "PTPN5", "PTPN6", "PTPN11", "PTPN21", "PTPN9", "PTPN2", "PTPN1")
    
    # Run UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
    
    # Generate plots
    PCA_plot <- PCAPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE)
    
    names(your_gene_list_PR) <- levels(seurat_obj)
    seurat_obj2 <- RenameIdents(seurat_obj, your_gene_list_PR)
    names(your_gene_list_PN) <- levels(seurat_obj)
    seurat_obj3 <- RenameIdents(seurat_obj, your_gene_list_PN)
    
    PR_cluster_annotate <- DimPlot(seurat_obj2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+labs(title = "PTP Receptors")
    PN_cluster_annotate <- DimPlot(seurat_obj3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+labs(title = "PTP non-Receptors")
    
    UMAP_plot <- UMAPPlot(seurat_obj,group.by = "seurat_clusters")
    
    PR_dotplot <- DotPlot(seurat_obj3, features = your_gene_list_PR, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(title = "PTP Receptors")
    PN_dotplot  <- DotPlot(seurat_obj, features = your_gene_list_PN, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(title = "PTP non-Receptors")
    
    PRfeat <- FeaturePlot(object = seurat_obj, features = your_gene_list_PR)
    PNfeat <- FeaturePlot(object = seurat_obj, features = your_gene_list_PN)
    
    PP_PN_cluster_annotate <- grid.arrange(PR_cluster_annotate, PN_cluster_annotate, ncol=2)
    
    PTFgeneexpression_dotplots <- grid.arrange(PR_dotplot, PN_dotplot, ncol=2)
    
    # Save the processed Seurat object for txt subset
    saveRDS(PCA_plot, "$pca_plot_file")
    saveRDS(PP_PN_cluster_annotate, "$pp_pn_cluster_annotate_file")
    saveRDS(UMAP_plot, "$umap_plot_file")
    saveRDS(PTFgeneexpression_dotplots, "$ptf_geneexpression_dotplots_file")
    saveRDS(PRfeat, "$pr_feat_file")
    saveRDS(PNfeat,  "$pn_feat_file")
    saveRDS(seurat_obj, "$final_seurat_obj_file")
    """
}

// Run the workflow
workflow {
    read_10x | preprocess_data
}
