#!/usr/local/bin/zsh
gosh -I. ica.scm | tee a.txt
gcsplit a.txt --prefix=t -b "%d.txt" '/^$/' '{*}'
rm t2.txt a.txt
mv *.txt text/
