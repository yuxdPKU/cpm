#!/bin/sh
set -e

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd "$srcdir"; \
aclocal -I "${OFFLINE_MAIN}/share"; \
libtoolize --force --copy; \
automake -a --add-missing --copy; \
autoconf)

"$srcdir/configure" "$@"
