#!/bin/bash

for sceneID in *
do


    cd $sceneID
    echo $sceneID

    img=`ls -1 | grep -v hdr`
    echo img $img
    hdr=`ls -1 | grep hdr`
    echo hdr $hdr

    new_img=$img"_MTLstack"
    new_hdr=`echo $hdr|sed 's/.hdr//'`
    new_hdr=$new_hdr"_MTLstack.hdr"

    echo img $img new_img $new_img
    echo hdr $hdr new_hdr $new_hdr
    mv       $img         $new_img
    mv       $hdr         $new_hdr

exit

    cd ..

done
