while read p; do
  echo "$p"
  wget http://www.proteinatlas.org/$p.json
done <prefixes.txt
