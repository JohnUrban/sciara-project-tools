##look at the distribution of overlap error rates in the unitig step to look for a peak in the distribution
##dump the overlap store yourself and compute the distribution using a command like on the unitigging or trimming stores
##This will dump overlaps for the first 1000 reads which should be a big enough sample to estimate the distribution. You'd then want to plot a histogram of the last column of the output.
~/software/canu/canu/Linux-amd64/bin/ovStoreDump -G g310-default-all.gkpStore -O g310-default-all.ovlStore -d -b 1 -e 1000 -coords -d5 -d3
