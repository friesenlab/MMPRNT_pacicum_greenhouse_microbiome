
#!/bin/bash

INDEX=1

 

for file_name in $(cat files.txt); do
  echo ${file_name}
  if ((${INDEX} % 2 == 0)); then

    echo "sample${INDEX},/data/friesen/MSU/alan_projects/LargeFiles_ForHPCC/20171213_16S-V4_PE/rawdata/${file_name},reverse"

  else

    echo "sample${INDEX},/data/friesen/MSU/alan_projects/LargeFiles_ForHPCC/20171213_16S-V4_PE/rawdata/${file_name},forward"

  fi

 

  INDEX=$((${INDEX} + 1))

done
