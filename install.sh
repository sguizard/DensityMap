#/bin/sh

mkdir /usr/local/share/DensityMap
cp ./* /usr/local/share/DensityMap
perl -i -pe 's/.\/color.txt/\/usr\/local\/share\/DensityMap\/color.txt/' /usr/local/share/DensityMap/DensityMap.pl
cd /usr/local/bin
ln -s /usr/local/share/DensityMap/DensityMap.pl

