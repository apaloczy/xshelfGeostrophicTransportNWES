for i in *.zip; do
  echo $i
  unzip $i
  rm $i
done
