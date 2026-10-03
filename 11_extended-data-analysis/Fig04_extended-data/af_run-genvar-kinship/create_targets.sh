awk '
{
  trait=$1
  for (mm_i in mms) {
    for (kin_i in kins) {
      for (idx_i in idxs) {
        print trait"_"mms[mm_i]"_"kins[kin_i]"_"idxs[idx_i]
      }
    }
  }
}
BEGIN {
  mms[1]="DMR"; mms[2]="SNP"
  kins[1]="SNP"; kins[2]="DMR"; kins[3]="genvar"
  idxs[1]="aBN"; idxs[2]="aIBS"
}
' traits > targets_all_DMR_SNP_all_kinships
