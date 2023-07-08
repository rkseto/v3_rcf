#!/bin/csh
rm -f tmpl.pid
echo R > tmpl.pid
./convert
rm -f tmpl.pid
echo $status > tmpl.pid
