# plot pie charts of superkingdom

# e.g.,
# R --slave --vanilla --args input.txt output samle_name < plot_pie.r

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

# trouble doing multiple pie charts in same picture ! 

args = commandArgs(TRUE);

inputfile=args[1];			# inputfile
outputfile=args[2];			# outputfile
samplename=args[3];			# name

if (file.info(inputfile)$size!=0)
{
	x=as.matrix(read.table(inputfile));
}

print(x);

png(filename=paste(outputfile,".png",sep=""),width=1024,height=640);

for (i in 1:length(x[,1])) 
{ 
	
	if (x[[i,1]]=="iteration") 
	{ 
		if (i!=1) 
		{
			print(superkingdom); 
			print(count); 
			pers=round(100*count[,]/sum(count)); 
			lbls <- paste(superkingdom, round(count), pers, sep=", "); 
			lbls <- paste(lbls,"%",sep=""); print(lbls);
			
			pie(count,labels = lbls, col=rainbow(length(lbls)),main=as.numeric(x[[i,2]])-1);
		}
		superkingdom=NULL; 
		count=NULL;
	} 
	else 
	{
		superkingdom=cbind(superkingdom,x[[i,1]]); 
		count=cbind(count,as.numeric(x[[i,2]]));
	} 
}

print(superkingdom); 
print(count); 
pers=round(100*count[,]/sum(count)); 
lbls <- paste(superkingdom, round(count), pers, sep=", "); 
lbls <- paste(lbls,"%",sep=""); print(lbls);

pie(count,labels = lbls, col=rainbow(length(lbls)),main="R1");

dev.off();
