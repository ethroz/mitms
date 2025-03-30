#!/bin/bash

sudo apt update
sudo apt install -y liblapack-dev gfortran

wget https://master.dl.sourceforge.net/project/lapackpp/lapackpp-2.5.4.tar.gz?viasf=1 -O lapackpp-2.5.4.tar.gz
tar -xf lapackpp-2.5.4.tar.gz
(
  cd lapackpp-2.5.4
  ./configure --prefix=/usr
  make
  sudo make install
)
rm -rf lapackpp-2.5.4*
