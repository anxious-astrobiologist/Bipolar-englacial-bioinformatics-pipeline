# This code was used to create figure 4 and was run using R v4.4.2 and R Studio v2024.12.0+467
library(dplyr)
library(readr)     # for parse_number()
library(ggplot2)
library(ggnewscale)

# 1) Ensure numeric columns are truly numeric (strip stray whitespace etc.)
# Input for this figure is a CSV file containing normalized metagenome and
# metatranscriptome counts of all marker genes as well as marker
# gene counts in each of the genome bins
Actic_Antarctic_normalized2_R <- Actic_Antarctic_normalized2_R %>%
  mutate(
    Metagenome.Counts  = parse_number(as.character(Metagenome.Counts)),
    Transcript.Counts  = parse_number(as.character(Transcript.Counts))
  )

# 2) Your ordering vectors (as in your code)
y_axis_order <- c(
  "NADH Oxidation (NuoF)", "Sulfide oxidation to sulfur (fccA)", 
  "Sulfide oxidation to sulfur (Sqr)", "Sulfide oxidtion to sulfite (dsrA)",
  "Thiosulfate oxidation (SoxB)", "Ammonia Oxidation (AmoA)", 
  "Nitrite Oxidation (NxrA)", "Iron (ll) oxidation (Cyc2)", 
  "Hydrogen Oxidation (hydA)", "Hydrogen Oxidation (hydB)", 
  "Anaerobic carbon monoxide oxidation (CooS)", 
  "Aerobic carbon monoxide oxidation (CoxL)", 
  "Methane Oxidation (PmoA)", "Methane Oxidation (mmoX)", 
  "Formate Oxidation (FdhA)", "Aerobic Respiration (high O2) (CoxA, CyoA)", 
  "Aerobic Respiration (low O2) (CydA, CcoN)", "Sulfate Reduction (DsrA)", 
  "Sulfate Reduction (AsrA)", "Fumerate Reduction (FrdA)", 
  "Nitrate Reduction (NarG)", "Nitrate Reduction (NapA)", 
  "Nitrite reduction to nitric oxide (NirS)", 
  "Nitrite reduction to nitric oxide (NirK)", 
  "Nitrite reduction to ammonium (NrfA)", "Nitric Oxide Reduction (NorB)", 
  "Nitrous Oxide Reduction (NosZ)", "Nitrogen Fixation (NifH)", 
  "Iron (lll) Reduction (MtrB)", "Reductive dehalogenation (Cbdb)", 
  "Methanogenesis (McrA)", "Hydrogen Production (hypF)", 
  "Oxygenic Photosynthesis (PsaA, PsbA)", 
  "Anoxygenic Phototrophy (PufM, PufL)", "Reductive TCA cycle (AclB)", 
  "Wood-Ljungdahl Pathway (AcsB)", "3-hydroxypropionate Cycle (Mcr)", 
  "Calvin-Benson-Bassham Cycle (RbcL)"
)
y_axis_reversed_order <- rev(y_axis_order)

x_axis_order <- c(
  "White Glacier", paste0("Bin ", 1:13, "_WG"),
  "Johnsons Glacier", paste0("Bin ", 1:17, "_JG")
)

# 3) Apply factor orders
Actic_Antarctic_normalized2_R$Metabolisms <- factor(
  Actic_Antarctic_normalized2_R$Metabolisms, levels = y_axis_reversed_order
)
Actic_Antarctic_normalized2_R$Sample <- factor(
  Actic_Antarctic_normalized2_R$Sample, levels = x_axis_order
)

# 4) Plot
heatmap <- ggplot(Actic_Antarctic_normalized2_R) + 
  geom_tile(aes(x = Sample, y = Metabolisms, fill = Metagenome.Counts)) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(1, 600),
    breaks = c(6, 60, 600),
    trans = "log",
    na.value = "white",
    name = "CPM"
  ) +
  labs(x = NULL, y = "Metabolisms") +
  theme_minimal() +
  
  # New fill scale for the points
  new_scale_fill() +
  geom_point(
    aes(x = Sample, y = Metabolisms, fill = Transcript.Counts),
    shape = 21, size = 5, colour = "black"
  ) +
  scale_fill_gradientn(
    colours = c("white", "red"),
    limits = c(1, 30000),
    breaks = c(3, 30, 300, 3000, 30000),
    trans = "log",
    na.value = "white",
    name = "TPM"
  ) +
  theme_minimal()
