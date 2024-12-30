export -f run_roast
run_roast() {
    roast -T={} -R 30 -M 20 -E Ath "((Aal Aar Aly Ane Aru Bca Bju Bna Bni Bol Bra Cen Chi Cru Csa Dni Ech Iin Lma Mpy Ovi Raq Rbr Rsa Sal Sar Spa Tar) Ath)" *.*.maf ref_species_multiway.maf > roast.sh
}
export -f run_roast

find . -name "*.maf" | parallel -j40 run_roast
