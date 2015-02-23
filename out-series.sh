#!/usr/local/bin/zsh
gosh -I. generate-random.scm > series.txt
gcsplit -b '%d.txt' series.txt '/^$/' '{*}'
cat xx0.txt | grep -v -e '^$' >  in.txt
cat xx1.txt | grep -v -e '^$' > out.txt
rm xx*.txt

