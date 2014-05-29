# Download some sample traces from FlyCircuit

# A set of 20 Kenyon cells used as a test data set in the nat package
# see ?kcs20 for details
kcs20_gene_names=c("FruMARCM-M001205_seg002", "GadMARCM-F000122_seg001", "GadMARCM-F000050_seg001", 
"GadMARCM-F000142_seg002", "FruMARCM-F000270_seg001", "FruMARCM-F001115_seg002", 
"FruMARCM-M001051_seg002", "GadMARCM-F000423_seg001", "ChaMARCM-F000586_seg002", 
"FruMARCM-M001339_seg001", "GadMARCM-F000476_seg001", "FruMARCM-F000085_seg001", 
"FruMARCM-F000706_seg001", "FruMARCM-M000842_seg002", "FruMARCM-F001494_seg002", 
"FruMARCM-F000188_seg001", "GadMARCM-F000071_seg001", "FruMARCM-M000115_seg001", 
"GadMARCM-F000442_seg002", "FruMARCM-F001929_seg001")

traceurls=file.path("http://flycircuit.tw","flycircuitImage","tracing",
  paste(names(kcs20),"_lineset.am.gz",sep=""))
  
if(!file.exists("tracings"))
  dir.create('tracings')

tracing_paths=file.path("tracings",basename(traceurls))
missing_traces=traceurls[!file.exists(tracing_paths)]

if(!length(missing_traces)) message("All tracings already downloaded!")
for (i in seq_along(missing_traces)){
  message("Downloading tracing ",i," out of ",length(missing_traces)," ...")
  download.file(url=missing_traces[i], destfile=file.path("tracings",basename(missing_traces[i])))
}
