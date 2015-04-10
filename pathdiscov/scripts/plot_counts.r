# plot the good blast hits, splits, and multiplicities, and on different tiers in the same graph

# e.g.,
# R --slave --vanilla --args input.txt output samle_name < plot_counts.r

args = commandArgs(TRUE);

inputfile=args[1];			# inputfile
outputfile=args[2];			# outputfile
samplename=args[3];			# name

# initialize vars
x=NULL;

if (file.info(inputfile)$size!=0)
{
	x=as.matrix(read.table(inputfile));
}

vertchunk=10;							# vertchunk defines how many pieces you cut the y-axis up into

y_min=0;
y_max=max(as.numeric(x[,2]));

png(filename=paste(outputfile,".png",sep=""),width=1024,height=640);

ydelta=(y_max-0)/vertchunk;

xrange=c(0,length(x[,2])+1);
yrange=c(0,y_max+ydelta);

# make a plot window
plot(xrange, yrange, type="n",xlab="filter iteration", ylab="# reads", xaxs = "i", yaxs = "i");
title(paste(samplename," - Iterative Read filtering",sep=""));

# add traces
lines(as.numeric(x[,2]), type="o", pch=16, lty=2, col="blue");

text(as.numeric(x[,2])+ydelta/2, x[,2], cex=0.8, col="red");

xvals=1:length(x[,2]);
yvals=rep(ydelta/4,length(x[,2]));

# text(xvals, yvals, x[,1], cex=1.0, col="black", srt=45);
text(xvals, yvals, x[,1], cex=0.8, col="black");

dev.off();