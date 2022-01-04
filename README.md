The article "Seriation using tree-penalized path length" by Aliyev & Zirbel, 2021, introduces a family of new, hybrid, seriation methods which utilize the strengths of two
well-known objectives, Path Length (TSP), and OLO (Optimal Leaf Order, a dendrogram seriation method by Bar-Joseph et al. 2001).
The main idea is to add to the distance matrix a penalty which measures the number and severity of deviations from the tree structure; this provides a balance between short path length and fidelity to the tree.
The new objective function is called tpPL, which stands for "tree-penalized path length".
Materials posted here include the datasets used in the study (44 datasets plus the R code that generated them), and example R files to give the user examples of using PL, OLO, and tpPL objective functions to produce orderings and plot the corresponding heat maps of the re-ordered distance matrix.
The user can choose the values for the following parameters:
b = (scaled) tree strength parameter;
link = type of clustering linkage of the hierarchical tree used to compute the tree penalty matrix;
rp = # of iterations for TSP and tpPL methods.

Visualizations of the datasets and the orderings produced for the paper are available at http://rna.bgsu.edu/publications/ordering/tppl/
