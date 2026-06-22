mkdir cellranger_out
mkdir cra_out
for LGnum in LG22 LG26 LG300 LG301 LG31 LG38; do
  cellranger-arc count --id=${LGnum}_HXB2 --reference=/projects1/references/CellRanger-GRCh38_HXB2/CellRanger-GRCh38_HXB2 --libraries=/projects1/newmulti/libraries/${LGnum}libraries.csv --localcores=16 --localmem=128 --create-bam=false
  mv ${LGnum} cellranger_out/
  ln -s ../cellranger_out/${LGnum}_HXB2/outs cra_out/${LGnum}
done
