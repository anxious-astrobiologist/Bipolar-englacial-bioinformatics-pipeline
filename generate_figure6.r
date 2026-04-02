#This code was used to generate Figure 6 using R v4.4.2 and R Studio v2024.12.0+467
library(ggplot2)
library(ggnewscale)

# Assuming you have a vector `y-axis_order` representing the desired order of levels in the y-axis
y_axis_order1 <- c("Cold shock, stress and HSP proteins", "Cold shock protein, CspA family", "Universal stress protein, UspA family", "Co-chaperonin GroES (HSP10)", "Chaperonin GroEL (HSP60 family)", "Molecular chaperone DnaK (HSP70)", "Molecular chaperone IbpA, HSP20 family", "Molecular chaperone, HSP90 family", "Molecular chaperone GrpE (heat shock protein)", "Ribosomal 50S subunit-recycling heat shock protein", "DNA replication and repair", "DNA gyrase/topoisomerase I, subunits A&B", "Recombinational DNA repair protein RecR", "RecA-superfamily ATPase, KaiC/GvpD/RAD55 family", "RecA/RadA recombinase", "DNA repair ATPase RecN", "Superfamily II DNA helicase RecQ", "Rad3-related DNA helicase", "Membrane and peptidoglycan alteration", "3-oxoacyl-(acyl-carrier-protein) synthase", "3-oxoacyl-[acyl-carrier-protein] synthase III (KASIII)", "3-oxoacyl-[acyl-carrier-protein] reductase", "Glycosyltransferase involved in cell wall biosynthesis", "3-hydroxyacyl-CoA dehydrogenase", "Fatty-acid desaturase", "D-alanyl-D-alanine carboxypeptidase", "Carotenoid biosynthesis", "Phytoene dehydrogenase-related protein", "Phytoene/squalene synthetase", "Polysaccharide capsule", "Capsular polysaccharide biosyn. protein EpsC", "Capsule polysaccharide export protein", "Capsular polysaccharide biosynthesis protein", "Exopolysaccharide biosynthesis protein", "Osmotic stress", "ABC proline/glycine betaine transport, ATPase component", "ABC proline/glycine betaine transport, permease", "Choline-glycine betaine transporter", "Osmoprotectant binding protein", "Choline dehydrogenase or related flavoprotein", "Trehalose-6-phosphate synthase", "Trehalose-6-phosphatase", "Maltooligosyltrehalose synthase", "Na+/proline symporter", "Na+/H+ antiporters", "Oxidative stress", "Catalase", "Peroxiredoxin", "Glutathione peroxidase", "Spermidine synthase", "Thioredoxin reductase", "Glyoxylase or related hydrolase, beta-lactamase superfam II", "Toxin/Antitoxin modules", "Ser/Thr kinase RdoA, MazF antagonist", "Antitoxin component of MazEF", "mRNA-interferase, toxin component of MazEF", "mRNA interferase YafQ, toxin component of YafQ-DinJ", "Antitoxin component of RelBE or YafQ-DinJ", "Translation and transcription factors", "tRNA-dihydrouridine synthase", "tRNA A37 threonylcarbamoyladenosine dehydratase", "Translation elongation factor EF-Tu, GTPase", "Translation elongation factor EF-G, GTPase", "Translation initiation factor IF-2, GTPase", "Translation initiation factor IF-3", "Transcription antitermination factor NusA", "Transcription termination factor NusB", "Superfamily II DNA and RNA helicase", "Superfamily II DNA or RNA helicase, SNF2", "Superfamily II RNA helicase")
# Reverse the order of `y_axis_order`
y_axis_reversed_order1 <- rev(y_axis_order1)
# Assuming you have a vector `x_axis_order` representing the desired order of levels in the x-axis
x_axis_order <- c("White Glacier", "Bin 1_WG", "Bin 2_WG", "Bin 3_WG", "Bin 4_WG", "Bin 5_WG", "Bin 6_WG", "Bin 7_WG", "Bin 8_WG", "Bin 9_WG", "Bin 10_WG", "Bin 11_WG", "Bin 12_WG", "Bin 13_WG", "Johnsons Glacier", "Bin 1_JG", "Bin 2_JG", "Bin 3_JG", "Bin 4_JG", "Bin 5_JG", "Bin 6_JG", "Bin 7_JG", "Bin 8_JG", "Bin 9_JG", "Bin 10_JG", "Bin 11_JG", "Bin 12_JG", "Bin 13_JG", "Bin 14_JG", "Bin 15_JG", "Bin 16_JG", "Bin 17_JG")  # Replace with your desired order

# Input file is a csv file containing normalized metagenome and metatranscriptome and bin counts of cold adapted marker genes

# Convert the `Metabolisms` column in your data frame to a factor with the desired order
Bipolar_cold_adapt_for_r$Gene <- factor(Bipolar_cold_adapt_for_r$Gene, levels = y_axis_reversed_order1)
# Convert the `Sample` column in your data frame to a factor with the desired order
Bipolar_cold_adapt_for_r$Sample <- factor(Bipolar_cold_adapt_for_r$Sample, levels = x_axis_order)

cold_adapt_heatmap <- ggplot(Bipolar_cold_adapt_for_r) + 
  geom_tile(aes(x = Sample, y = Gene, fill=Metagenome.Counts)) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(1, 6000),
    breaks = c(6, 60, 600, 6000),
    trans = "log",
    na.value = 'white',
    name = "CPM"
  ) +
  labs(x = element_blank(), y = 'Cold Response Genes', fill = 'CPM') +
  theme_minimal() +
  # Command for a new scale
  new_scale_fill() +
  geom_point(aes(x = Sample, y = Gene, fill = Transcript.Counts), shape = 21, size = 5, colour = "black") +
  scale_fill_gradientn(
    colours = c("white", "red"),
    limits = c(1, 10500),
    breaks = c(10, 100, 1000, 10000),
    trans = 'log',
    na.value = 'white') +
  #scale_fill_viridis_c(
  #  option = "plasma",
  #  limits = c(1, 30000),
  #  breaks = c(3, 30, 300, 3000, 30000),
  #  trans = 'log',
  #  na.value = 'white'
  #) +
  labs(x = element_blank(), y = 'Cold Response Genes', fill = 'TPM') +
  theme_minimal()
