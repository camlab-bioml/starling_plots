rm(list=ls())
library(CATALYST)
library(SingleCellExperiment)
library(reticulate)
np <- import("numpy")
library(anndata)

cohort <- 'basel'
seg_type <- 'dc'

fns <- list.files(paste0("/home/campbell/yulee/project/st/", cohort, "/mat/", seg_type))
file_ids <- sort(as.integer(do.call(c, strsplit(fns, ".npy"))))

df <- NULL
for ( file_id in file_ids ) {
  tmp <- np$load(paste0("/home/campbell/yulee/project/st/", cohort, "/mat/", seg_type, "/", file_id, ".npy"))
  df <- rbind(df, tmp)
}

if (cohort == 'tonsil') {
  
  df <- df[df[,1] > 3,]
  df[,1] <- df[,1] - 3
  
  info <- c("SMA(Y89Di)",           "Y89Di",     "Y89",
            "E_Cadherin(In113Di)",  "In113Di",    "In113",    
            "Cytokeratin(In115Di)", "In115Di",    "In115",
            "HLA_DR(Pr141Di)",     "Pr141Di",    "Pr141",
            "Vimentin(Nd142Di)",   "Nd142Di",    "Nd142",    
            "CD28(Nd144Di)",        "Nd144Di",    "Nd144",    
            "CD15(Nd145Di)",        "Nd145Di",    "Nd145",    
            "CD45RA(Nd146Di)",      "Nd146Di",    "Nd146",    
            "CD66b(Sm147Di)",       "Sm147Di",    "Sm147",    
            "CD20(Sm149Di)",        "Sm149Di",    "Sm149",    
            "CD68(Nd150Di)",        "Nd150Di",    "Nd150",    
            "CD4(Eu151Di)",        "Eu151Di",    "Eu151",    
            "CD8(Sm152Di)",         "Sm152Di",    "Sm152",    
            "CD11c(Sm154Di)",       "Sm154Di",    "Sm154",    
            "CD45RO(Dy162Di)",      "Dy162Di",    "Dy162",    
            "CD3(Dy163Di)",         "Dy163Di",    "Dy163",    
            "IFNg(Tm169Di)",        "Tm169Di",    "Tm169",    
            "TCF1(Er170Di)",        "Er170Di",    "Er170",    
            "CD14(Yb172Di)",        "Yb172Di",    "Yb172",    
            "CD56(Yb173Di)",        "Yb173Di",    "Yb173",    
            "PD1(Yb176Di)",         "Yb176Di",    "Yb176", 
            "DNA1(Ir191)",          "Ir191Di",    "Ir191",
            "DNA2(Ir193)",          "Ir193Di",    "Ir193",          
            "CD45(Pt194Di)",        "Pt194Di",    "Pt194",    
            "PNAd(Pt195Di)",        "Pt195Di",    "Pt195",    
            "CD31(Pt196Di)",        "Pt196Di",    "Pt196")
} else if (cohort == 'basel') {
  info <- c(
    "Argon(Ar80)", "Ar80Di", "Ar80",
    "RT2(Ru96)", "Ru96Di", "Ru96",
    "RT3(Ru98)", "Ru98Di", "Ru98",
    "RT4(Ru99)", "Ru99Di", "Ru99",
    "RT5(Ru100)", "Ru100Di", "Ru100",
    "RT6(Ru101)", "Ru101Di", "Ru101",
    "RT7(Ru102)", "Ru102Di", "Ru102",
    "RT8(Ru104)", "Ru104Di", "Ru104",
    "RT9(Ru10)", "Ru10Di", "Ru10",
    "totHH3(In113)", "In113Di", "In113",
    "H3K27me3(La139)", "La139Di", "La139",
    "CK5(Pr141)", "Pr141Di", "Pr141",        
    "Fibronectin(Nd142)", "Nd142Di", "Nd142",
    "CK19(Nd143)", "Nd143Di", "Nd143", 
    "CK8n18(Nd144)", "Nd144Di", "Nd144",
    "Twist(Nd145)", "Nd145Di", "Nd145",
    "CD68(Nd146)", "Nd146Di", "Nd146",
    "CK14(Sm147)", "Sm147Di", "Sm147",
    "SMA(Nd148)", "Nd148Di", "Nd148",
    "Vimentin(Sm149)", "Sm149Di", "Sm149",
    "cMyc(Nd150)", "Nd150Di", "Nd150",
    "HER2(Eu151)", "Eu151Di", "Eu151",
    "CD3(Sm152)", "Sm152Di", "Sm152",
    "pHH3(Eu153)", "Eu153Di", "Eu153",
    "ERK(Sm154)", "Sm154Di", "Sm154",
    "Slug(Gd155)", "Gd155Di", "Gd155",
    "RabbitIgGHL(Gd156)", "Gd156Di", "Gd156",
    "PR1(Gd158)", "Gd158Di", "Gd158",
    "PR2(Gd158)", "Gd158Di", "Gd158",
    "p53(Tb159)", "Tb159Di", "Tb159",
    "CD44(Gd160)", "Gd160Di", "Gd160",
    "EpCAM(Dy161)", "Dy161Di", "Dy161",      
    "CD45(Dy162)", "Dy162Di", "Dy162",
    "GATA3(Dy163)", "Dy163Di", "Dy163",
    "CD20(Dy164)", "Dy164Di", "Dy164",
    "BCatenin(Ho165)", "Ho165Di", "Ho165",    
    "CAIX(Er166)", "Er166Di", "Er166",
    "ECadherin(Er167)", "Er167Di", "Er167",
    "Ki67(Er168)", "Er168Di", "Er168",
    "EGFR(Tm169)", "Tm169Di", "Tm169",
    "S6(Er170)", "Er170Di", "Er170",
    "Sox9(Yb171)", "Yb171Di", "Yb171",
    "vWF(Yb172)", "Yb172Di", "Yb172",
    "CD31(Yb172)", "Yb172Di", "Yb172",
    "mTOR(Yb173)", "Yb173Di", "Yb173",
    "CK7(Yb174)", "Yb174Di", "Yb174",
    "pCK(Lu175)", "Lu175Di", "lu175",
    "kE(Lu175)", "Lu175Di", "Lu175",
    "cPARP(Yb176)", "Yb176Di", "Yb176",
    "cCaspase3(Yb176)", "Yb176Di", "Yb176",
    "DNA1(Ir191)", "Ir191Di", "Ir191",
    "DNA2(Ir193)", "Ir193Di", "Ir193")
  
  #Lu175,pan Cytokeratin,AE1,142,0.5 ug/mL,0.352112676,1,20,46
  #Lu175,Keratin Epithelial, clone AE3,200,0.5 ug/mL,0.25,0,,46
  #Yb176,cleaved PARP,F21-852,100,8 ug/mL,8,1,,47
  #Yb176,Cleaved Caspase3,C92-605,150,8 ug/mL,5.333333333,0,,47
  
} else if (cohort == 'meta') {
  info <- c(
    "Argon(Ar80)", "Ar80Di", "Ar80",
    "totHH3(In113)", "In113Di", "In113",
    "Xe126(Xe126)", "Xe126Di", "Xe126",
    "I127(I127)", "I127Di", "I127",
    "Xe131(Xe131)", "Xe131Di", "Xe131",
    "Xe134(Xe134)", "Xe134Di", "Xe134",
    "H3K27me3(La139)", "La139Di", "La139",
    "Ce140(Ce140)", "Ce140Di", "Ce140",
    "CK5(Pr141)", "Pr141Di", "Pr141",
    "Fibronectin(Nd142)", "Nd142Di", "Nd142",
    "CK19(Nd143)", "Nd143Di", "Nd143",
    "CK8n18(Nd144)", "Nd144Di", "Nd144",
    "Twist(Nd145)", "Nd145Di", "Nd145",
    "CD68(Nd146)", "Nd146Di", "Nd146",
    "CK14(Sm147)", "Sm147Di", "Sm147",
    "SMA(Nd148)", "Nd148Di", "Nd148",
    "Vimentin(Sm149)", "Sm149Di", "Sm149",
    "cMyc(Nd150)", "Nd150Di", "Nd150",
    "HER2(Eu151)", "Eu151Di", "Eu151",
    "CD3(Sm152)", "Sm152Di", "Sm152",
    "pHH3(Eu153)", "Eu153Di", "Eu153",
    "ERK(Sm154)", "Sm154Di", "Sm154",
    "Slug(Gd155)", "Gd155Di", "Gd155",
    "ER(Gd156)", "Gd156Di", "Gd156",          
    "PR(Gd158)", "Gd158Di", "Gd158",
    "p53(Td159)", "Td159Di", "Td159",
    "CD44(Gd160)", "Gd160Di", "Gd160",
    "EpCAM(Dy161)", "Dy161Di", "Dy161",
    "CD45(Dy162)", "Dy162Di", "Dy162",
    "GATA3(Dy163)", "Dy163Di", "Dy163",
    "CD20(Dy164)", "Dy164Di", "Dy164",
    "BCatenin(Ho165)", "Ho165Di", "Ho165",
    "CAIX(Er166)", "Er166Di", "Er166",
    "ECadherin(Er167)", "Er167Di", "Er167",
    "Ki67(Er168)", "Er168Di", "Er168",
    "EGFR(Tm169)", "Tm169Di", "Tm169",
    "pS6(Er170)", "Er170Di", "Er170",
    "Sox9(Yb171)", "Yb171Di", "Yb171",
    "vWFCD31(Yb172)", "Yb172Di", "Yb172",
    "mTOR(Yb173)", "Yb173Di", "Yb173",
    "CK7(Yb174)", "Yb174Di", "Yb174",
    "pCK(Lu175)", "Lu175Di", "Lu175",
    "cPARP(Yb176)", "Yb176Di", "Yb176",
    "DNA1(Ir191)", "Ir191Di", "Ir191",
    "DNA2(Ir193)", "Ir193Di", "Ir193",
    "Hg202(Hg202)", "Hg202Di", "Hg202",
    "Pb204(Pb204)", "Pb204Di", "Pb204",
    "Pb206(Pb206)", "Pb206Di", "Pb206",
    "Pb207(Pb207)", "Pb207Di", "Pb207",
    "Pb208(Pb208)", "Pb208Di", "Pb208")
}

info <- as.data.frame(matrix(info, ncol = 3, byrow = TRUE))
rownames(info) <- info[,1]; info <- info[,-1]
colnames(info) <- c("channel_name", "marker_name")

if (seg_type == 'dc') {
  
  sm <- read.csv("/home/campbell/yulee/project/st/spillover.csv")
  rownames(sm) <- sm[,1]; sm <- sm[,-1]
  
  X <- df[,9:ncol(df)]
  sce1 <- SingleCellExperiment(list(exprs=t(X)))
  rowData(sce1) <- info
  sce1 <- compCytof(sce1, sm, assay = "exprs", method = "nnls", transform = FALSE, cofactor = 1, overwrite = FALSE)
  
  out <- cbind(df[,1:8], t(assays(sce1)$compcounts))
} else {
  out <- cbind(df[,1:7], 0, df[,8:ncol(df)])
}

colnames(out) <- c('sample_id', 'cell_id', 'x', 'y', 'area', 'area_convex', 'has_neighbor', 'nuclei_count', do.call(rbind, strsplit(rownames(info), "\\("))[,1])

ad <- AnnData(
  X = out[,9:ncol(out)],
  obs = as.data.frame(out[,1:8]),
  var = data.frame(marker = colnames(out)[9:ncol(out)]))

write.csv(as.data.frame(out), file = paste0("/home/campbell/yulee/project/st/", cohort, "/exp_mat/", seg_type, "/exp_mat.csv"))
write_h5ad(ad, paste0("/home/campbell/yulee/project/st/", cohort, "/exp_mat/", seg_type, "/exp_mat.h5ad"))