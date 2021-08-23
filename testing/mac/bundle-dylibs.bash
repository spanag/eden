#!/bin/bash

set -e


# from https://stackoverflow.com/a/21188136
get_abs_filename() {
  # $1 : relative filename
  if [ -d "$(dirname "$1")" ]; then
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  fi
}

# Set up parameters
if ( [ -z "$1" ] )||( [ -z "$2" ] ); then
echo "Usage: $0 <executable to get dylibs for> <path tree to copy linked dylibs from>"
exit 1
fi

FIX_EXE=$(get_abs_filename "$1")
REDIST_PATH=$(get_abs_filename "$2/")
EXE_PATH=`dirname $FIX_EXE`
# echo redist $REDIST_PATH 
# echo exe    $EXE_PATH 

if ( [ -z "$DYLIB_SUBDIR"] ); then
DYLIB_SUBDIR=dylibs
fi
# Note: $DYLIB_SUBDIR is a relative path to the exe. Get the relative path between exe and dylibs more generally,
# with this subroutine: https://unix.stackexchange.com/a/85068 

DYLIB_DIR="$EXE_PATH/$DYLIB_SUBDIR"
mkdir -p "$DYLIB_DIR"


# otool -L "$FIX_EXE" 
# echo

# crop all characters before and after... with .*before and after
item_list_str=`otool -L "$FIX_EXE" | \
sed -n "s,.*$REDIST_PATH/\(.*dylib\).*,\1,p"`

# Copy the dylibs to the new subdir, and collect the linker flags
edit_flags=
IFS=$'\n';
for item in $item_list_str
do
	# unset IFS
	
	item_dir=`dirname $item`
	item_base=`basename $item`
	
	# echo $item $item_dir $item_base
	
	# Place all dependent dylibs in a hierarchy, may result in big paths
	# mkdir -p "$DYLIB_DIR/$item_dir"
	# cp -f "$REDIST_PATH/$item" "$DYLIB_DIR/$item"
	# edit_flags+=" -change $REDIST_PATH/$item @executable_path/$DYLIB_SUBDIR/$item"
	
	# Alternatively:
	# Place all so's in one level, hope there are no name conflicts
	cp -f "$REDIST_PATH/$item" "$DYLIB_DIR/$item_base"
	edit_flags+=" -change $REDIST_PATH/$item @executable_path/$DYLIB_SUBDIR/$item_base"
	
	# echo
	
	# IFS=$'\n';
done
unset IFS

if ( [ -n "$edit_flags" ] ); then
echo "Fixing references:"
echo install_name_tool $edit_flags "$FIX_EXE"
install_name_tool $edit_flags "$FIX_EXE"
fi


