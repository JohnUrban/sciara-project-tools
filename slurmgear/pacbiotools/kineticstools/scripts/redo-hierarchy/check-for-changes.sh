tail contig*/kin*csv > t1
sleep 1
tail contig*/kin*csv > t2
diff t1 t2 | grep "^>"
