#/bin/sh

mkdir /usr/local/share/DensityMap
cp -r ./* /usr/local/share/DensityMap
perl -i -pe 's/.\/colours.txt/\/usr\/local\/share\/DensityMap\/colours.txt/' /usr/local/share/DensityMap/DensityMap.pl
cd /usr/local/bin
ln -s /usr/local/share/DensityMap/DensityMap.pl

