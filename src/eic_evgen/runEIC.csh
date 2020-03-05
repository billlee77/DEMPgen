#!/bin/csh


rm -f ./eic
g++ -o eic eic.cxx `root-config --cflags --glibs`

set j = 0
set pol = 1
set type = 2
# pol = 1 ( Polarization Up )  pol = 2 ( Polarization Down )
while ( $j < 1 )

   << EOF ./eic
$pol
$type
1000000
$j
EOF
      
   @ j++
end

