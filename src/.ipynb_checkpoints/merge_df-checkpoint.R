library(reshape2)
# FUNCTION TO CREATE nxp MATRIX
# WITH n=number_patients
# AND p include gene level features
merge_clinical_mutation <- function(dd_clinical, dd_maf, row_field="TARGET_NAME", col_field="GENE", binary=TRUE) {
    # INPUT:
    # dd_clinical clinical dataframe nxq
    # dd_maf      mutation maf dataframe mx
    # binary: if true then binary answer per gene, if false then number of alterations in that gene
    # OUTPUT:
    # list with:
    # - ddmerge concatenates clinical and mutations 
    # - ddmut mutations only
    library(reshape2)
    # reshape
    tmp.binary  = dcast(data=dd_maf, as.formula(paste0(row_field,"~",col_field)), value.var=row_field, fun.aggregate=length)
    rownames(tmp.binary) = tmp.binary[,1]
    tmp.binary = tmp.binary[,-1]
    # binary answer
    if (binary) {
        tmp.binary[tmp.binary > 0] = 1
    }
    # merge
    ddm = dd_clinical
    ddm[,colnames(tmp.binary)] = 0
    im = match(ddm$LEUKID, rownames(tmp.binary))
    ddm[!is.na(im),colnames(tmp.binary)] = tmp.binary[im[!is.na(im)],]
    # mutations only
    ddmut = ddm[,((ncol(ddm)-ncol(tmp.binary)+1):ncol(ddm))]
    rownames(ddmut) = ddm$LEUKID

    return(list(ddmerge=ddm, ddmut=ddmut))
}
