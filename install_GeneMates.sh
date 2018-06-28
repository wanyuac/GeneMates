#!/bin/bash
# Copyright 2017-2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 2 July 2017; lastest edition: 17 Janurary 2018

display_usage() {
    echo "Installing GeneMates on a Linux machine.
    This script uses at most two positional parameters.
    Usage: (bash) ./install_GeneMates.sh [package] [install to which directory]
    Examples:
        chmod u+x install_GeneMates.sh
        ./install_GeneMates.sh GeneMates_0.1.0.tar.gz ~/R_lib
        ./install_GeneMates.sh ~/R_lib  # install the default package (specified in this script) to a target directory
        ./install_GeneMates.sh  # install the default package under the default directory of R libraries specified by .libPaths()[1].
    "
}

# Set default values for user variables
package="GeneMates.tar.gz"
exit_code=0

# Read arguments and make an installation of the package
if [ $# -eq 0 ]; then  # only $0
    if [ -f $package ]; then
        echo "Since neither the package nor the library has been specified, the default package ${package} will be installed to the GeneMates folder under the default directory of R libraries."
        R CMD INSTALL $package
    else
        echo "Error: the package ${package} does not exist. Installation failed."
        display_useage
        exit_code=126  # The command is not executable.
	fi
elif [ $# -eq 1 ]; then  # $0 and $1
    lib_dir=$1
    if [ ! -d $lib_dir ]; then
        echo "Error: the directory ${lib_dir} does not exist. Installation failed."
        display_usage
        exit_code=126
    elif [ -f $package ]; then
        echo "Since no package has been specified, the default package ${package} will be installed to the GeneMates folder under the directory ${lib_dir}."
        R CMD INSTALL --library=$lib_dir $package
    else
        echo "The package ${package} does not exist. Installation failed."
        display_usage
        exit_code=126  # The command is not executable.
	fi
else  # at least three parameters
    package=$1
    lib_dir=$2
    if [ ! -d $lib_dir ]; then
        echo "Error: the directory ${lib_dir} does not exist. Installation failed."
        display_usage
        exit_code=126
    elif [ -f $package ]; then
        echo "The package ${package} will be installed to the GeneMates folder under the directory ${lib_dir}"
        R CMD INSTALL --library=$lib_dir $package
    else
        echo "Error: the package ${package} does not exist. Installation failed."
        display_usage
        exit_code=126  # The command is not executable.
	fi
fi

exit $exit_code
