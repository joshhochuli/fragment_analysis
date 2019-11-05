#!/bin/bash

for i in *
do 
  move=False
  for j in $i/*
  do
    count=$(wc -l $j | awk -F " " '{print $1}')
    if [ $count == 0 ]
    then
      move=True
      echo $j
    fi
  done
  if [ $move == True ]
  then
    echo $i
    mv $i dirs_with_issues/empty_files/.
  fi
done
