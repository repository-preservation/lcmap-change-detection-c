#!/bin/bash

########################################################################
#
# tileLandsat.sh
#
# Parent script to "tile" landsat scenes into files of stacked bands
# in ENVI format, in BIP order, sequenced as required by the CCDC.
#
# Assumes inputs are in a directory and ESPA-packaged .tar.gz files.
# Output is written to sub-directories under a parent directory,
# one per scene ID.
#
# 20151019 bdavis
# Original development, plagiarizing  from fireMapping.sh.  Input
# scene IDs should be arguments, somehow.
# 
########################################################################


########################################################################
#
# Set up environment.  This should be valid for any of the SLURM nodes
# in the EROS YetiJr environment.
#
########################################################################

export PATH=.:/alt/local/run:/alt/local/bin:/usr/local/local.host/bin:/usr/lib64/qt-3.3/bin:/usr/local/local.host/sbin:/sbin:/bin:/usr/sbin:/usr/bin:/root/bin:/home/bdavis/bin
export LD_LIBRARY_PATH=/alt/local/lib:/usr/local/local.host/ssl/lib:/usr/local/local.host/lib64:/usr/local/local.host/lib:/lib64:/usr/lib64:/usr/local/local.host/mysql/lib
export PYTHONPATH=/alt/local/bin/:/alt/local/lib/python2.7/site-packages


########################################################################
#
# Parse arguments.  So far just in and out paths.
#
########################################################################

date
echo ""

if [ $# -ne 3 ]; then
    echo "Usage: tileLandsat.sh input-path output-path exe"
    echo "Where input-path is the directory containing ESPA-created .tar.gz files,"
    echo "and output-path is the parent directory in which to write results in sub-directories."
    echo "and exe is the type of executable to prepare for, MatLAB or C."
    echo ""
    echo "For Example: tileLandsat.sh /data/bdavis/WA/2005-2015/Albers_grid07 /data/bdavis/WA/2005-2015/Tiles_grid07 C"
    echo ""
    exit
fi


#in=/data/bdavis/WA/2005-2015/Albers_grid07
#out=/data/bdavis/WA/2005-2015/Tiles_grid07

in=$1
out=$2
exe=$3
originalDir=`pwd`

echo Reading from $in, writing to $out
echo ""
#echo in $in
#echo out $out
#echo originalDir $originalDir

########################################################################
#
# It seems one could just read from /hsm, if it were mounted, but it's
# not.
#
# Also, the 20% clear requirement is outdated thinking, because WA data
# is ESPA rectangles (Area of Interest), with much more temporal density
# (and fill), but also the change for any one x/y pixel location to be
# clear/snow/water, even though the total scene/ESPA-AOI is less than
# 20%.  Therefore, process without checking for clear threshold and
# skipping, but for the C verion only. Retain the 20% restriction for
# the MatLAB version, because that is what is expected.
#
########################################################################

########################################################################
#
# ls -1 LE7*.gz|sort -n --key=1.10,1.16
# will give all path/row file names in a gridNN directory, sorted by
# date, if that is required.  Zhe says no, his software sorts
# internally.checking with Song.
# There are 174 46/27 LE7 in grid07, n of them are clear enough.
#
# This list of scenes needs to be parameterized.  Piping an ls,
# cat-ing a text file, etc. Limted to the first 24 path 46 row 27
# LE7 files for initial testing purposes.
#
########################################################################

#for sceneID in LE70460272005010 LE70460272005042 LE70460272005058 LE70460272005074 LE70460272005106 LE70460272005154 LE70460272005170 LE70460272005186 LE70460272005202 LE70460272005218 LE70460272005234 LE70460272005250 LE70460272005266 LE70460272005282 LE70460272005298 LE70460272005314 LE70460272005330 LE70460272005346 LE70460272006045 LE70460272006061 LE70460272006077 LE70460272006109 LE70460272006125 LE70460272006141

cd $in

#just zz scenes
#for sceneID in `ls -1 LE7046027*.gz|sort -n --key=1.10,1.16`
#production
for sceneID in `ls -1 *.gz|sort -n --key=1.10,1.16`
#testing
#for sceneID in `ls -1 LT4*.gz|sort -n --key=1.10,1.16`
#for sceneID in `ls -1 LT5*.gz|sort -n --key=1.10,1.16`
#for sceneID in `ls -1 LE7*.gz|sort -n --key=1.10,1.16`
#for sceneID in `ls -1 LC8*.gz|sort -n --key=1.10,1.16`
#for sceneID in LT50450281996227 # test of zero valid pixels
#for sceneID in LT50450281996243 # test of 100 pct cloud cover
#for sceneID in LT50450281995112 # test of valid scene

    do

    ####################################################################
    #
    # Set up the names to use, and unpackage the scene.
    #
    ####################################################################

    sceneID=`echo $sceneID|cut -c 1-16`
    pkg=`ls -1 $in/*.gz|grep $sceneID`
    echo sceneID $sceneID
    echo pkg $pkg
    echo ""

    # if the directory and hdr already exist, assume this is done, skip.
    if [ -d $out/$sceneID ]; then

        hdr=`ls -1 $out/$sceneID/*.hdr|wc|awk '{print $1}'`
        echo hdr
        ls $out/$sceneID/*hdr

        if [ $hdr -ne 0 ]; then
            echo skipping hdr
            ls $out/$sceneID
            echo ""
            continue
        fi
    fi

    mkdir -p $out/$sceneID
    tar -C $out/$sceneID -zxvf $pkg
    cd $out/$sceneID
    name=`ls -1 *cfmask.tif|cut -c 1-21`
    echo name $name



    ####################################################################
    #
    # At some point, we may need to extract the sensor, to determine
    # sensor-specific processing of which bands, etc.
    #
    ####################################################################

    sensor=`echo $sceneID|cut -c 1-3`


    ####################################################################
    #
    # If the scene is not at least 20% clear, skip.  Clear is defined
    # as clear + water / total.  Call a perl tool modified from 
    # fireMapping.sh which calculated cloud cover, because we need
    # clear.  FYI, "clear" is a reserved word, hence "cleer".
    # Attempting to do floating point math in bash was causing some
    # divide by zero errors because of rounding very small values
    # (less than 0.01) to integer values.  cloudCover.pl also returns
    # zero for scenes whose values are all fill, so zero valid pixels.
    # (yes, I've found one.)
    #
    ####################################################################

    ccargs=`gdalinfo -hist *cfmask.tif|tail -2|head -1|awk '{print $1 " " $2 " " $3 " " $4 " " $5}'`
    echo ccargs $ccargs
    pctClear=`cloudCover.pl $ccargs`
    intPctClear=`echo $pctClear | cut -d '.' -f 1`

    echo Percent clear pixels in total valid pixels: $intPctClear
    echo ""

    if [ $intPctClear -lt 20 ] && [ $exe == "MatLAB" ]; then
        echo Percent clear pixels less than 20: $intPctClear
        echo ""
        #rm -f *.tif *.xml *.txt
        cd ..
        rm -f -r $sceneID
        cd $originalDir
        continue
    fi


    ####################################################################
    #
    # For the Song C version of ccdc, individual envi files for each
    # band are required.  For the Zhe MatLAB version, stacked band envi
    # files in BIP format are required.
    #
    ####################################################################

    if [ $exe == "C" ]; then

        if [ $sensor == "LC8" ]; then

            gdal_translate -of ENVI $name"_sr_band2.tif"   $name"_sr_band2.img"
            gdal_translate -of ENVI $name"_sr_band3.tif"   $name"_sr_band3.img"
            gdal_translate -of ENVI $name"_sr_band4.tif"   $name"_sr_band4.img"
            gdal_translate -of ENVI $name"_sr_band5.tif"   $name"_sr_band5.img"
            gdal_translate -of ENVI $name"_sr_band6.tif"   $name"_sr_band6.img"
            gdal_translate -of ENVI $name"_sr_band7.tif"   $name"_sr_band7.img"
            gdal_translate -of ENVI $name"_toa_band10.tif" $name"_toa_band10.img"
            gdal_translate -of ENVI $name"_cfmask.tif"     $name"_cfmask.img"

            mv *.img *sr_band2.hdr $out/.

        else

            gdal_translate -of ENVI $name"_sr_band1.tif"  $name"_sr_band1.img"
            gdal_translate -of ENVI $name"_sr_band2.tif"  $name"_sr_band2.img"
            gdal_translate -of ENVI $name"_sr_band3.tif"  $name"_sr_band3.img"
            gdal_translate -of ENVI $name"_sr_band4.tif"  $name"_sr_band4.img"
            gdal_translate -of ENVI $name"_sr_band5.tif"  $name"_sr_band5.img"
            gdal_translate -of ENVI $name"_sr_band7.tif"  $name"_sr_band7.img"
            gdal_translate -of ENVI $name"_toa_band6.tif" $name"_toa_band6.img"
            gdal_translate -of ENVI $name"_cfmask.tif"    $name"_cfmask.img"

            mv *.img *sr_band1.hdr $out/.

        fi

        cd ..
        rm -r $sceneID
        cd $originalDir

    else


        ################################################################
        #
        # Convert the cfmask from byte to 16bit.  All "bands" in the
        # merged tif need to be of the same data type.
        #
        ################################################################

        cfmask=`ls -1 *_cfmask.tif`
        cfmask16=`echo $cfmask | sed 's/cfmask/cfmask16/'`
        echo cfmask $cfmask cfmask16 $cfmask16
        echo ""
        echo "Converting cfmask to 16-bit."
        echo ""
        gdal_translate -ot UInt16 $cfmask  $cfmask16
        echo ""


        ################################################################
        #
        # Merge all the bands together, in the specified order, in ENVI
        # format.
        #
        ################################################################

        echo "Merging bands."
        echo ""

        if [ $sensor == "LC8" ]; then

            gdal_merge.py -o $name"_stack.img"  \
                          -separate             \
                          -of ENVI              \
                          $name"_sr_band2.tif"  \
                          $name"_sr_band3.tif"  \
                          $name"_sr_band4.tif"  \
                          $name"_sr_band5.tif"  \
                          $name"_sr_band6.tif"  \
                          $name"_sr_band7.tif"  \
                          $name"_toa_band10.tif" \
                          $name"_cfmask16.tif"

        else

            gdal_merge.py -o $name"_stack.img"  \
                          -separate             \
                          -of ENVI              \
                          $name"_sr_band1.tif"  \
                          $name"_sr_band2.tif"  \
                          $name"_sr_band3.tif"  \
                          $name"_sr_band4.tif"  \
                          $name"_sr_band5.tif"  \
                          $name"_sr_band7.tif"  \
                          $name"_toa_band6.tif" \
                          $name"_cfmask16.tif"

        fi

        echo ""


        ################################################################
        #
        # Create a jpg just during testing if a sanity check is
        # required.
        #
        ################################################################

#        echo "Creating JPEG."
#        echo ""
#        gdal_translate -of JPEG          \
#                       -b 6 -b 4 -b 3    \
#                       -ot Byte          \
#                       -scale            \
#                       $name"_stack.img" \
#                       $name.jpg
#        echo ""


        ################################################################
        #
        # Translate the envi file to BIP.  This could be combined with
        # the initial merge if we do not need a bsq jpg.  Eventually.
        #
        ################################################################

        echo "Converting to BIP."
        echo ""
        gdal_translate -co "INTERLEAVE=BIP" \
                       -of ENVI             \
                       $name"_stack.img"    \
                       $name"_MTLstack.img"
        echo ""


        ################################################################
        #
        # Clean up after yourself.  Leaves only the required files and
        # reduces confusion.  (was going to say "eliminates".)
        #
        ################################################################

        mv $name"_MTLstack.img" $name"_MTLstack"
        rm *.tif *.xml *_stack.hdr *_stack.img
        cd ..

        echo ""

    fi

    done

echo "Processing complete."
echo ""
date

exit

