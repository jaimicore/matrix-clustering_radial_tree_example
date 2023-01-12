# ========================== #
# Run RSAT matrix-clustering #
# ========================== #

# We assume you already have installed RSAT in your system or alternatively you can run it in the webserver: http://rsat.sb-roscoff.fr/matrix-clustering_form.cgi

# matrix-clustering -v 2                                            \
# -matrix JASPAR2022_CORE_nematodes JASPAR2022_CORE_nematodes.tf tf \
# -title 'JASPAR 2022 nematodes CORE'                               \
# -hclust_method average                                            \
# -calc sum                                                         \
# -metric_build_tree Ncor                                           \
# -lth w 5 -lth cor 0.6 -lth Ncor 0.4                               \
# -label_in_tree name                                               \
# -return json                                                      \
# -quick                                                            \
# -radial_tree_only                                                 \
# -o results/JASPAR2022_CORE_nematodes/JASPAR2022_CORE_nematodes


# =============== #
# Add annotations #
# =============== #

## Adds the annotation +  color ring + background
Rscript annotate_matrix-clustering.R                                                     \
	--annotation data/Extra/JASPAR_2022_nematodes_CORE_annotations_radial_tree.tsv   \
	--input      data/JASPAR_2022_CORE_nematodes_matrix-clustering
