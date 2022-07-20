

#   scale_fill_manual(values = friendly_cols)
#   scale_color_manual(values = friendly_cols)

library(ggplot2)

theme_multipanel =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(size=0.1),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		legend.position = "bottom",
		strip.background = element_blank(),
		axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
		panel.spacing.x=unit(0.1, "lines"),
		axis.text.x = element_text(size=6),
		axis.text.y = element_text(size=6),
		strip.text.x = element_text(size = 7),
		strip.text.y = element_text(size = 7),
		
		# legend
		legend.key.size = unit(5, 'mm'),
		legend.title = element_text(size=7),
		legend.text = element_text(size=6),
		
		# Avoid text clipping for facets. Currently not merged remotes::install_github("tidyverse/ggplot2#4223")
		#strip.clip = "off",
		
		# Title
		plot.title = element_text(size=7),
		
		axis.line.x = element_line(size=0.2),
		axis.line.y = element_line(size=0.2),
		axis.ticks.x = element_line(size=0.2),
		axis.ticks.y = element_line(size=0.2)
	)

# Patchwork
# +  plot_layout(guides = 'collect' ) + plot_annotation(tag_levels = c('A')) & theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'))

save(theme_multipanel, file="data/theme_multipanel.rda", compress = "xz")
