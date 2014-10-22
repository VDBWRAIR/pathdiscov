# plot pie charts of superkingdom

# WORST SCRIPT EVER ! REWRITE INTO SOMETHING DECENT !

# e.g.,
# R --slave --vanilla --args input.txt output sample_name < plot_phylo_percents.r

# input file should be of the form, e.g.:

# iteration	1
# Viruses	1479.26
# Eukaryota	124.88
# Unannotated	325.84
# Bacteria	437.96
# iteration	2
# Viruses	1739.26
# Eukaryota	391.4
# Unannotated	18.24
# iteration	3
# Viruses	88.51
# Eukaryota	11.5
# Bacteria	0.1
# Unannotated	258.9

args = commandArgs(TRUE);

inputfile=args[1];			# inputfile
outputfile=args[2];			# outputfile
samplename=args[3];			# name

if (file.info(inputfile)$size!=0)
{
	x=as.matrix(read.table(inputfile));
}

print(x);

y_min=0;
y_max=100;

x_min=1;
x_max=0;

line_vir=NULL;
line_euk=NULL;
line_bac=NULL;
line_un=NULL;

png(filename=paste(outputfile,".png",sep=""),width=1024,height=640);

for (i in 1:length(x[,1])) 
{ 
	
	if (x[[i,1]]=="iteration") 
	{ 
		if (i!=1) 
		{
#			print(superkingdom); 
#			print(count); 
			pers=round(100*count[,]/sum(count)); 
			lbls <- paste(superkingdom, pers);
			print(lbls); 			 
			
			# HORRIBLE STYLE but use it for now:
			for (j in 1:length(superkingdom)) 
			{
				if (superkingdom[[j]]=="Viruses") {line_vir[x_max]=pers[[j]]}
				if (superkingdom[[j]]=="Eukaryota") {line_euk[x_max]=pers[[j]]}
				if (superkingdom[[j]]=="Bacteria") {line_bac[x_max]=pers[[j]]}
				if (superkingdom[[j]]=="Unannotated") {line_un[x_max]=pers[[j]]}
			}
		}
		superkingdom=NULL; 
		count=NULL;
		x_max=x_max+1;
	} 
	else 
	{
		superkingdom=cbind(superkingdom,x[[i,1]]); 
		count=cbind(count,as.numeric(x[[i,2]]));
	} 
}

#print(superkingdom); 
#print(count); 
pers=round(100*count[,]/sum(count)); 
lbls <- paste(superkingdom, pers);
print(lbls);  

# HORRIBLE STYLE but use it for now:
for (j in 1:length(superkingdom)) 
{
	if (superkingdom[[j]]=="Viruses") {line_vir[x_max]=pers[[j]]}
	if (superkingdom[[j]]=="Eukaryota") {line_euk[x_max]=pers[[j]]}
	if (superkingdom[[j]]=="Bacteria") {line_bac[x_max]=pers[[j]]}
	if (superkingdom[[j]]=="Unannotated") {line_un[x_max]=pers[[j]]}
}

xrange=c(x_min,x_max);
yrange=c(y_min,y_max);

# print(line_vir);
# print(line_euk);
# print(line_bac);
# print(line_un);

# make a plot window
plot(xrange, yrange, type="n",xlab="iteration", ylab="%");
title(paste(samplename," - Iterative Blast",sep=""));

# add traces
lines(line_vir, type="o", pch=16, lty=2, col="blue");
lines(line_euk, type="o", pch=16, lty=2, col="red");
lines(line_bac, type="o", pch=16, lty=2, col="green");
lines(line_un, type="o", pch=16, lty=2, col="purple");

text(cbind(1,1,1,1), cbind(2+line_vir[1], 2+line_euk[1], 2+line_bac[1], 2+line_un[1]), cbind("Viruses", "Eukaryota", "Bacteria", "Unannotated"), cex=0.8, col=cbind("blue","red","green","purple"));

dev.off();
