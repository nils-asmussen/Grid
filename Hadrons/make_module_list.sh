#!/usr/bin/env bash

echo 'modules_cc =\' > modules.inc
find Modules -name '*.cc' -type f -print | sort \
   | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
echo 'modules_hpp =\' >> modules.inc
find Modules -name '*.hpp' -type f -print | sort \
   | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
for f in `find Modules -name '*.hpp'`; do
   echo "#include <Hadrons/${f}>"
done | sort > Modules.hpp
