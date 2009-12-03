#!/bin/bash
BOOST_PREFIX=/home/gav/Projects/Software/usr

echo "This will create a statically built version of qctool, called qctool-static, in the current directory."
echo "Please set $BOOST_PREFIX in the script according to your machine."

g++ -O3 -o qctool-static src/*.cpp genfile/src/*.cpp statfile/src/*.cpp $BOOST_PREFIX/lib/libboost_iostreams-gcc43-mt-1_39.a $BOOST_PREFIX/lib/libboost_filesystem-gcc43-mt-1_39.a /usr/lib/libz.a -static-libgcc -static -I include/ -I genfile/include/ -I statfile/include/ -I $BOOST_PREFIX/include/boost-1_39/ /usr/lib/libz.a -I build/release/statfile/ -DHAVE_ZLIB=1 -DHAVE_BOOST=1 -DHAVE_BOOST_FILESYSTEM=1 -I build/release/genfile/include/ apps/qctool.cpp $BOOST_PREFIX/lib/libboost_system-gcc43-mt-1_39.a


