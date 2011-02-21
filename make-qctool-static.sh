#!/bin/bash
#
# This script creates a fully static build of qctool.
# It is intended for use on Linux machines; on Mac OS, one should only really link to static
# boost libraries which can be done with some juggling of the waf configure options.
# For example, link to the static iostreams, filesystem, and system libs from a new directory boost-static,
# and run waf configure --static --boost-libs=/Users/gav/Projects/Software/qctool/boost-static
# The result can be checked with otool build/release/qctool

BOOST_PREFIX=/home/gav/Projects/Software/usr
TARGET=qctool-static-linux-x86-64
echo "This will create a statically built version of qctool, called $TARGET, in the current directory."
echo "Please set \$BOOST_PREFIX in the script according to your machine."
echo "I have it as \"$BOOST_PREFIX\"."

g++ -O3 -o $TARGET src/*.cpp string_utils/src/*.cpp genfile/src/*.cpp statfile/src/*.cpp appcontext/src/*.cpp fputils/src/*.cpp worker/src/*.cpp $BOOST_PREFIX/lib/libboost_iostreams-gcc43-mt-1_44.a $BOOST_PREFIX/lib/libboost_filesystem-gcc43-mt-1_44.a /usr/lib/libz.a -static-libgcc -static -I include/ -I string_utils/include/ -I genfile/include/ -I statfile/include/ -I appcontext/include/ -I fputils/include/ -I worker/include/ -I $BOOST_PREFIX/include/boost-1_44/ /usr/lib/libz.a -I build/release/statfile/ -DHAVE_ZLIB=1 -DHAVE_BOOST=1 -DHAVE_BOOST_FILESYSTEM=1 -I build/release/genfile/include/ apps/qctool.cpp $BOOST_PREFIX/lib/libboost_system-gcc43-mt-1_44.a $BOOST_PREFIX/lib/libboost_thread-gcc43-mt-1_44.a -lpthread

mv $TARGET /tmp/qctool
chmod ugo+rwx /tmp/qctool
gzip /tmp/qctool
mv /tmp/qctool.gz ./$TARGET.gz


