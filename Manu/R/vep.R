# VEP INPUT

load("/data/xwang/adsp3/sampling.rdt")
mcmc_id = gsub("_.*", "", rownames(mcmc)) # 9726

vep.input = read.table("/data/xwang/adsp3/plink/wgs.bim")
vep.input = vep.input[vep.input$V2 %in% mcmc_id, ]

vep.input$V7 = paste(vep.input$V5, vep.input$V6, sep = "/")
vep.input = vep.input[c("V1", "V4", "V4", "V7")]
vep.input$V8 = "+"

file = "~/Dropbox/GitHub/wgs2/Manu/R/vep_input.txt"
write.table(vep.input, file = file, quote = F, sep = "\t", row.names = F)

file = "~/Dropbox/GitHub/wgs2/Manu/R/vep_output.txt"
vep = read.table(file = file, stringsAsFactors = F, header = T)

save(vep, file = "~/Dropbox/GitHub/wgs2/Manu/R/vep.rdt")
