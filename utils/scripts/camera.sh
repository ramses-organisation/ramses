for i in `\ls -d output*`; do cat ${i}/halo_0*.txt0* | sort -r -nk 2 | head -n 1; done

