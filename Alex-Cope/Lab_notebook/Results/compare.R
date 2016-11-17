plotComp <- function(together,cat_1,cat_2,title,main_1,main_2,ylab_1,ylab_2)
{
	pdf(title)
	sel <- read.table(together,sep=",",header=TRUE)
	sel_c1 <- read.table(cat_1,sep=",",header=TRUE)
	sel_c2 <- read.table(cat_2,sep=",",header=TRUE)
	plot(sel[,3],sel_c1[,3],main=main_1,xlab="Together",ylab=ylab_1)
	ribModel:::upper.panel.plot(sel[,3],sel_c1[,3],sd.x=sel[,4:5],sd.y=sel_c1[,4:5])
	plot(sel[,3],sel_c2[,3],main=main_2,xlab="Together",ylab=ylab_2)
	ribModel:::upper.panel.plot(sel[,3],sel_c2[,3],sd.x=sel[,4:5],sd.y=sel_c2[,4:5])
	dev.off()
}

# plotComp("../Nosp_w_mp/selection.csv","../Nosp_w_mp_categories/selection_nosp.csv",
# 	"../Nosp_w_mp_categories/selection_mp.csv",
# 	"Nosp_vs_mp.pdf",
# 	"Together vs Genes w/ No SigPep",
# 	"Together vs Just Mature Peptides",
# 	"Genes w/ No SigPep",
# 	"Mature Peptides")
# 
# 
# plotComp("../Nosp_w_sp/selection.csv","../Nosp_w_sp_categories/selection_nosp.csv",
# 	"../Nosp_w_sp_categories/selection_sp.csv",
# 	"Nosp_vs_sp.pdf",
# 	"Together vs Genes w/ No SigPep",
# 	"Together vs Just Signal Peptides",
# 	"Genes w/ No SigPep",
# 	"Signal Peptides")

plotComp("../Mp_w_sp/Run_Sim/selection_mp_w_sp.csv","../Mp_w_sp_cat/First_run_sim/selection_mp.csv",
	"../Mp_w_sp_cat/First_run_sim/selection_sp.csv",
	"Mp_vs_sp_sim.pdf",
	"Together vs Mature Peptides",
	"Together vs Signal Peptides",
	"Mature Peptides",
	"Signal Peptides")