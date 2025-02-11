

## load
# ps.glom.norm.top
# from D6336 and D6300
# use plotfun
# and merge plots
# save --> Fig X.



### TO DO:
# source the codes in the chunks into different .R files!
scales::show_col(pal.man)

cowplot::plot_grid(
  scales::show_col(brewer.pal(10, 'Set1')),
  scales::show_col(brewer.pal(10, 'Accent')),
  scales::show_col(brewer.pal(10, 'Set2')),
  scales::show_col(brewer.pal(10, 'Dark2')),
  ncol = 1)



colors   <- setdiff(as.character(ps.glom.norm.top.merged@tax_table[,t]), 'unassigned')

unknowns <- grep("unknown_", colors, value = T)
unknowns <- as.character(grey.colors(length(unknowns)))
names(unknowns) <- grep("unknown_", colors, value = T)
unknowns <- c(unknowns, unassigned = 'black')

colors <- grep("unassigned", colors, value = T, invert = T)

colors <- grep("unknown_", colors, value = T, invert = T)

my.pal <- pal.man[1:length(colors)]

names(my.pal) <- colors

my.pal <- c(my.pal, unknowns)

scales::show_col(my.pal)
scales::show_col(pal.man)
#factor()

