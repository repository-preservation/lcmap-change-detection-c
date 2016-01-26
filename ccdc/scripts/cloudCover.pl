#!/usr/bin/perl

# ######################################################################
#
# Name: cloudCover.pl
#
# Author: bdavis
#
# Date: 20151029
#
# Description:
# Copied from cc.pl.  Receives as arguments the output from gdal -hist
# the first 5 values of which are number of clear, water, cloud shadow,
# snow, and cloud, respectively.  It returns the percentage of clear
# pixels as determined by clear + water / total.  cc.pl returns percent
# cloud cover.  This is called from bash scripts because floating
# point math is difficult in bash, and rounding errors were causing
# divide by zero when calculating the total pixels as a percentage of
# possible pixels in a grid.
#
# As it turns out, we have come across 1 ESPA scene which was all fill
# values, so zero valid pixels, resulting in divide by zero, hence the
# check for total not equal to zero.
# 
# Multipy result times by 100 because tileLandsat.sh is expecting
# percents as integer values.
#
# ######################################################################

$clear       = $ARGV[0];
$water       = $ARGV[1];
$cloudShadow = $ARGV[2];
$snow        = $ARGV[3];
$cloud       = $ARGV[4];
##print "clear         $clear\n";
##print "water         $water\n";
##print "cloudShadow   $cloudShadow\n";
##print "snow          $snow\n";
##print "cloud         $cloud\n";

$total = $clear + $water + $cloudShadow + $snow + $cloud;
#$cc = (($cloudShadow + $cloud) / $total ) * 100;
if ($total != 0)
    {
    $pctClear = (($clear + $water) / $total ) * 100;
    }
else
    {
    $pctClear = 0;
    }

## debug
# test of 1 for total gave a value of 4e-08 for pctTotal,
# which still was interpreted as a valid non-zero value
# and did not give a divide by zero error.
#$pctTotal = $total / (5000 * 5000);
#print "pctTotal $pctTotal\n";
#print "total $total\n";
#print "cc $cc\n";
## debug
#print "$cc\n";
print "$pctClear\n";

exit;
