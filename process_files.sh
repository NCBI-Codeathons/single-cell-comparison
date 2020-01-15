# Output the */* portion matching and then run the cut command on the files themselves
for FPATH in $(find /data/counts/*/*.counts); do
  # cut off .counts ending and return the group cell combination
  PROJECT_CELL=$(echo $FPATH | cut -d '/' -f 4,5 | cut -d '.' -f 1)
  # Get specific GROUP
  PROJECT=$(echo $PROJECT_CELL |  cut -d '/' -f 1)
  # Get CELL
  CELL=$(echo $PROJECT_CELL |  cut -d '/' -f 2)
  mkdir -p /data/gene_count_length_files/$PROJECT
  touch /data/gene_count_length_files/$PROJECT/$CELL.tsv
  cut -f 1,6,7 /data/counts/$PROJECT/$CELL.counts > /data/gene_count_length_files/$PROJECT/$CELL.tsv
done
