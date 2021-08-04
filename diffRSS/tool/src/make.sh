#!/bin/bash

javac -d ./ -cp ./lib/* bed/*.java count/*.java exonseq/*.java genome/*.java ./main/*.java mapping/*.java sim/*.java split/*.java
jar -cvmf MANIFEST.MF ../../diffRSS.jar *