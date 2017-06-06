xset <- xcmsSet(step=0.02,snthresh=3,mzdiff = 0.05)
grp <- group(xset,bw=10,mzwid = 0.05)
pc.tmp = ""
an <- annotate(grp, cor_eic_th=0, cor_exp_th=0,na.ok=TRUE)
write.csv(an$annotated,'data0210mzwid05diff05.csv')

