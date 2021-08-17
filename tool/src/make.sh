#!/bin/bash

javac -d ./ -cp lib/htsjdk-2.10.1.jar:lib/commons-math3-3.6.1.jar:lib/com.github.haifengl.smile-math-2.6.0.jar bed/*.java count/*.java exonseq/*.java genome/*.java ./main/*.java mapping/*.java sim/*.java split/*.java
jar -cvmf MANIFEST.MF ../../diffRSS.jar *
