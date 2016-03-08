#!/bin/bash

source ./version.in

version_day=`date +"%d"`
version_month=`date +"%m"`
version_year=`date +"%Y"`
version_yy=`date +"%y"`

#version_major=1
#version_minor=0
#version_build=9
#version_revis=2
#version_cont=0

if [ "$#" -eq 0 ] ; then
  version_cont=$(($version_cont+1))
else
  if [ "$1" == "-r" ] ; then
    version_cont=0
    version_revis=$(($version_revis+1))
  fi
fi

version_hist=$(($version_hist+1))
cat > include/version.h << EOF

#ifndef VERSION_H
#define VERSION_H

namespace AutoVersion{
	
	//Date Version Types
	static const char DATE[] = "${version_day}";
	static const char MONTH[] = "${version_month}";
	static const char YEAR[] = "${version_year}";
	static const char UBUNTU_VERSION_STYLE[] = "${version_yy}.${version_month}";
	
	//Software Status
	static const char STATUS[] = "Alpha";
	static const char STATUS_SHORT[] = "a";
	
	//Standard Version Type
	static const long MAJOR = ${version_major};
	static const long MINOR = ${version_minor};
	static const long BUILD = ${version_build};
	static const long REVISION = ${version_revis};
	
	//Miscellaneous Version Types
	static const long BUILDS_COUNT = ${version_cont};
	#define RC_FILEVERSION ${version_major},${version_minor},${version_build},${version_revis}
	#define RC_FILEVERSION_STRING "${version_major}, ${version_minor}, ${version_build}, ${version_revis}\\${version_cont}"
	static const char FULLVERSION_STRING [] = "${version_major}.${version_minor}.${version_build}.${version_revis}";
	
	//These values are to keep track of your versioning state, don't modify them.
	static const long BUILD_HISTORY = ${version_hist};
	

}
#endif //VERSION_H

EOF

cat > version.in << EOF
# current version
version_day=$version_day
version_month=$version_month
version_year=$version_year
version_yy=$version_yy
version_major=$version_major
version_minor=$version_minor
version_build=$version_build
version_revis=$version_revis
version_cont=$version_cont
version_hist=$version_hist
EOF
