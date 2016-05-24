#/bin/sh

mkdir /usr/local/share/DensityMap
cp -r ./* /usr/local/share/DensityMap
cd /usr/local/bin
ln -s /usr/local/share/DensityMap/DensityMap.pl

