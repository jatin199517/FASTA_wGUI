## Starting claudeR MCP
# Install
if (!require("devtools")) install.packages("devtools")
## devtools::install_github("IMNMV/ClaudeR")
## For vscode only
## devtools::install_github("nx10/httpgd")

# Set up your AI tool (For RStudio)
##library(ClaudeR)
##install_clauder()          # For Claude Desktop / Cursor
##install_cli(tools = "claude")  # For Claude Code CLI

## Add new message to the conversation
##claudeAddMessage("Hello, I'm analyzing flow cytometry data for a CAR T cell project")
# Start the server in RStudio
#claudeAddin()
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
library(purrr)

# Reading metadata file
listOfLibrary <- c("this.path","dplyr","tidyr","ggthemes","readxl","writexl","rstudioapi","ggplot2",
                   "flowCore","flowStats","ggcyto","openCyto","flowWorkspace","patchwork",
                   "ggpubr","flowMeans","grDevices","future","furrr","devtools",
                   "BiocManager","concaveman","stringr","rstudioapi","lubridate","knitr","kableExtra","httpgd")

# Installing BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# GitHub fallbacks for packages not on CRAN/Bioconductor (pkg name -> "user/repo").
# Add an entry here for any GitHub-only package added to listOfLibrary above.
githubRepos <- c(
  httpgd  = "nx10/httpgd",
  ClaudeR = "IMNMV/ClaudeR"
)

# Install/load each package, cascading across all 3 repos: CRAN -> Bioconductor
# -> GitHub. requireNamespace() is re-checked between tiers because
# install.packages() only warns (does not error) when a package is unavailable
# on CRAN, so relying on the error handler alone would skip the fallbacks.
walk(listOfLibrary, ~{
  if (!requireNamespace(.x, quietly = TRUE)) {
    message("Trying CRAN for: ", .x)
    tryCatch(install.packages(.x), error = function(e) NULL)
  }
  if (!requireNamespace(.x, quietly = TRUE)) {
    message("CRAN failed for: ", .x, ". Trying Bioconductor...")
    tryCatch(BiocManager::install(.x, ask = FALSE, update = FALSE),
             error = function(e) NULL)
  }
  if (!requireNamespace(.x, quietly = TRUE)) {
    repo <- githubRepos[[.x]]
    if (is.null(repo))
      stop("Package '", .x, "' not found on CRAN or Bioconductor, and no GitHub ",
           "repo mapping exists. Add it to githubRepos as ", .x, " = \"user/repo\".")
    message("Bioconductor failed for: ", .x, ". Trying GitHub: ", repo)
    devtools::install_github(repo)
  }
  library(.x, character.only = TRUE)
})
setwd(tryCatch(
  this.dir(),
  error = function(e) dirname(sub("^file://", "", utils::URLdecode(this.path::this.path(original = TRUE))))
))
message("[FASTA] Working directory: ", getwd())
FunctionFile<-grepv("ExtFunctions",list.files(pattern=".*\\.R$"))
listOfFunctions<-c(FunctionFile)
# Bootstrap: locate and source ExtFunctions.R, then let its .get_script_dir()
# helper set the working directory authoritatively. A minimal locate step is
# unavoidable here -- we must find ExtFunctions.R before we can source the
# helper that normalizes the path. In a VSCode R terminal the wd is already the
# workspace folder; in RStudio we fall back to the active document's folder.

walk(listOfFunctions, source)          # defines .get_script_dir() + pipeline functions

# Writing basic metadata file and allowing the user add further data.

metadata<-data.frame(path=list.files(pattern="\\.fcs$", recursive=TRUE, full.names = TRUE)) %>% 
  mutate(fcs_filename=(basename(path)),
         Year=str_extract(path,"/(\\d{4})/", group=1),
         SampleType=str_extract(fcs_filename, "PBMC"),
         DateOfAcquisition=sapply(path, function(f) d<-keyword(read.FCS(f, transformation=FALSE, which.lines=1,truncate_max_range = F))[["$DATE"]]),
         Compensation_type=NA, Downsampling=NA, Downsampling_n=NA) %>% mutate(StainType=str_extract(fcs_filename,"Full|FM.*(?=_Samples)|Unstained"),
                                                                              label=paste0(fcs_filename)) %>% 
  
  mutate(GatingMarker=str_extract(StainType,"FM[^_]+_(.+)", group=1)) %>% 
  dplyr::filter(!str_detect(StainType, "Unstained"))

cs<-load_cytoset_from_fcs(metadata$path)
markers<-data.frame(dims=colnames(cs), markers=sapply(colnames(cs), function(x){
  index<-which(names(markernames(cs)) %in% x)
  if(length(index)==0) return(NA_character_) else return(unname(markernames(cs)[index]))
}),
TransformType=NA, TransformArgs=NA)

# Retrieving updated instruction file from user.
instructions_filename<-list.files(pattern="instructions+\\.xlsx")
listOfSheets<-excel_sheets(instructions_filename)
finalInstructions<-setNames(map(listOfSheets,~read_excel(instructions_filename, sheet=.x)),listOfSheets)
## Changing markernames
updatedMarkernamesIndex<-which(!is.na(finalInstructions$marker$markers))
updatedMarkernames<-finalInstructions$marker$markers[updatedMarkernamesIndex]
names(updatedMarkernames)<-finalInstructions$marker$dims[updatedMarkernamesIndex]
markernames(cs)<-updatedMarkernames

##pData assignment
merged<-merge(pData(cs)[, "name", drop=FALSE], metadata, by.x="name", by.y="fcs_filename") %>% relocate(name)
rownames(merged)<-merged$name
pData(cs)<-merged

## Compensation
cs_comp <- apply_compensation(cs, metadata = metadata)

## Transformation
cs_comp_trans <- apply_transformation(cs_comp, markers = finalInstructions$markers)

## Gating
gs <- apply_gating2(cs_comp_trans, template = finalInstructions$gating_template, mc.cores = 4)

# Auto layouts
plot_list <- apply_plotting(gs, finalInstructions$gating_template)

# See what's available
names(plot_list)                        # sample names
names(plot_list[[2]])                   # layout names for first sample

# Access by index
plot_list[[7]][["layout_1"]]
# View all layouts for one sample
##for (layout in plot_list[[2]]) print(layout)

# Panel QC [Inactivated for Github]
##PanelMasterInventoryLocation<-"/Users/jsharma@lyell.com/Library/CloudStorage/OneDrive-LyellImmunopharma/Research - Bioinformatics & Tumor-Tcell Profiling/Flow Cytometry/TIER 1/InventoryManagement/260519_Master_Panel_Inventory.xlsx"
##PanelIteration<-"260527_LYL273_[JS]"
##pqc <- panel_QC(gs, PanelMasterInventoryLocation, PanelIteration)

# Gating QC
qc <- run_gating_qc(gs, finalInstructions$gating_template)
tbls <- gating_qc_tables(qc)

## Stats
pDataForJoining<-pData(gs) %>% mutate(sample=rownames(pData(gs)))
stats<-gs_pop_get_stats(gs,type="Percent") %>%
  left_join(gs_pop_get_stats(gs,type="Count")) %>% 
  left_join(pDataForJoining) %>% 
  group_by(sample) %>% 
  mutate(pop=str_remove(pop,"/notDebris/")) %>% 
  mutate(`QC: Count > 40`=ifelse(count>40,"Yes","No")) %>% 
  mutate(percent=100*percent)
## Summary plots
date<-paste0(str_extract(now(),"(\\d{4})-(\\d{2})-(\\d{2})",group=c(1,2,3)), collapse="")
##BID<-unique(stats$BenchlingID)
subtitle<-paste0("ELN ID"," | PBMC immunophenotyping assay")
SummaryPlot1<-stats %>% 
  dplyr::filter(str_detect(pop,"/CD4.CD8.")) %>% 
  ggplot(aes(x = `StainType`, y = percent, fill=`QC: Count > 40`)) +
  facet_wrap(~pop, scale="free_x") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ##stat_summary(aes(group = StudyID), geom = "line", size = 0.5, linetype = "solid") + # Connect the means
  geom_point(alpha = 0.7, size = 3, stroke = 0.8, shape = 21) +  # Enhanced settings
  scale_y_continuous(breaks=seq(0,100,20),expand=expansion(mult=c(0.1,0.1)),minor_breaks = seq(0,100,5)) +
  labs(
    title = paste0("Data Summary: 1"),
    subtitle=paste0(subtitle," | Populations"),
    y="% of Parent",
    x="FDP_lot") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size=13),
        axis.text.y = element_text(size=13),
        ##panel.grid.minor = element_line(color = "#4F5259", linewidth = 0.1),
        axis.title=element_text(size=15, face="bold"),
        legend.title=element_text(size=12),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        panel.border = element_rect(color = "#4F5259", fill = NA, linewidth = 0.8), 
        strip.text = element_text(size = 7, face = "bold", color="white"),
        plot.title = element_text(size = 18, face = "bold"),
        strip.background = element_rect(fill = "maroon"),
        panel.spacing = unit(0.4, "lines"),
        plot.subtitle = element_text(size = 13, color = "gray50", face="italic"),
        plot.margin = margin(t = 5, r = 20, b = 30, l = 5, unit = "pt")
  ) + coord_cartesian(clip = "off")
SummaryPlot1

## Export stats
date<-paste0(str_extract(now(),"(\\d{4})-(\\d{2})-(\\d{2})",group=c(1,2,3)), collapse="")
##write_xlsx(list("Stats"=stats), paste0(date,"_",BID,"_LYL273_cPARP.xlsx"))






