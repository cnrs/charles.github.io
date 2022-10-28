Identification and Validation of Aging-Related Genes in Idiopathic Pulmonary Fibrosis  
https://www.frontiersin.org/articles/10.3389/fgene.2022.780010/full  


Identification and Validation of Aging-Related Genes in Alzheimer’s Disease  
https://www.frontiersin.org/articles/10.3389/fnins.2022.905722/full  


Human Aging Genomic Resources (https://genomics.senescence.info/)
```
wget https://genomics.senescence.info/genes/human_genes.zip
```

MSigDB gene sets (http://www.gsea-msigdb.org/gsea/msigdb/index.jsp) The gene sets of "c5.go.bp.v7.5.1.symbols" was downloaded from the MSigDB database
```
wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/c5.go.bp.v2022.1.Hs.symbols.gmt
```

GO:0007568 aging (http://amigo.geneontology.org/amigo/term/GO:0007568)  
GO:0090398 cellular senescence (http://amigo.geneontology.org/amigo/term/GO:0090398)
```
grep AGING c5.go.bp.v2022.1.Hs.symbols.gmt
GOBP_AGING      http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_AGING KRTAP4-9        CDKN2B  CARM1   GNA13   CTSC    GJB6    SEC63   CNP     CNR1    COL4A2  COMP       TPRA1   ADM     SLC32A1 CRYAB   KRT25   ADRA1A  CCN2    CTNNA1  CYP1A1  ADRB3   DAG1    ACE     DDC     NQO1    DNMT3A  ABAT    AGT     EDN1    EDNRA   EDNRB   EIF2S1     ELAVL4  ENDOG   EPO     ERCC1   ERCC2   FGF2    FOXG1   SCAP    FOXO3   DNMBP   FOS     ALPL    NR5A1   SIN3A   FBXO4   AMFR    GFRA1   AMH     GJB2    SERP1   GCLC       GCLM    HTRA2   GNA11   GNA12   GPX1    GPX4    GRB2    GSK3A   APAF1   APEX1   HTR2A   APOD    IGF1R   IGFBP1  IGFBP2  IGFBP5  IL10    IL15    INHBA   INHBB   IRAK1      JUND    KCNMB1  ARG1    KRT14   KRT16   KRT33B  KRT83   HELT    LOXL2   SMAD4   MBP     FOXO4   MME     MMP2    MPO     ASS1    MT-ATP6 MT-CO1  NUDT1   MT-ND4  NFE2L2     NPY     NPY5R   ATP2B1  NTRK1   RNF165  OGG1    ATP5F1A P2RY1   PAX2    PAX5    GLRX2   PDGFRB  SERPINF1        ATP8A2  PENK    PITX3   MBD3    POLB    AVPR1A  AVPR1B     PPP3CA  RETN    PTGS2   BAK1    HAMP    RPN2    RPS6KB1 BGLAP   CX3CL1  SRR     SLC1A2  SLC6A3  SLC12A2 SOD1    SOD2    SREBF1  TACR3   TFRC    TGFB3   TGFBR2  THTSPO     TIMP1   SERPING1        C1QA    TNFRSF1B        KRTAP4-8        TYMS    UCP2    UCP3    VCAM1   WRN     CALCA   CTC1    CANX    CASP2   CASP9   PPP1R9B CAT     KRTAP4-5   KRTAP4-3        KMO     TP63    BECN1   HYAL2   EIF2B5  MBD2    CLDN1   RGN     AURKB   LONP1   KL      CD68    KCNE2
GOBP_MULTICELLULAR_ORGANISM_AGING       http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_MULTICELLULAR_ORGANISM_AGING  GNA13   SEC63   COMP    CRYAB   ADRA1A  DDEDN1     EDNRA   ERCC1   ALPL    NR5A1   SERP1   GNA11   GNA12   IGFBP1  INHBA   INHBB   HELT    RNF165  AVPR1A  AVPR1B  SLC1A2  TH      WRN     TP63    HYAL2


grep SENESCENCE c5.go.bp.v2022.1.Hs.symbols.gmt
GOBP_CELLULAR_SENESCENCE        http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_CELLULAR_SENESCENCE   AKT3    MIR543  CDK2    CDK6    CDKN1A  ZMPSTE24        CDKN1B     CDKN2A  CDKN2B  CITED2  KAT5    PLK2    NEK6    ZNF277  CGAS    MAPK14  VASH1   PLA2R1  SMC5    SIRT1   MORC3   NUP62   ABL1    ULK3    RSL1D1  FBXO5   MAGEA2B NSMCE2     H2AX    HLA-G   HMGA1   HRAS    ID2     IGF1R   ING2    KIR2DL4 ARG2    LMNA    ARNTL   MIR10A  MIR146A MIR17   MIR188  MIR217  MIR22   MIR34A  MAGEA2  MAP3K3  MAP3K5     MIF     MNT     ATM     NPM1    YBX1    OPA1    PAWR    ABI3    FZR1    WNT16   SIRT6   PML     PRMT6   PRELP   SLC30A10        PRKCD   MAPK8   MAPK11  MAPK9   MAPK10     MAP2K1  MAP2K3  MAP2K6  MAP2K7  B2M     ZMIZ1   PTEN    MIR20B  RBL1    BCL6    MAP2K4  BMPR1A  SPI1    SRF     BRCA2   NEK4    TBX2    TBX3    MIR590  TERC    TERF2      TERT    TP53    TWIST1  WNT1    WRN     SMC6    KAT6A   ZKSCAN3 HMGA2   CALR    YPEL3   ECRG4   MAPKAPK5        TP63    PNPT1   DNAJA3  EEF1E1  NUAK1
GOBP_REPLICATIVE_SENESCENCE     http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_REPLICATIVE_SENESCENCE        CDKN1A  CDKN2A  CHEK1   CHEK2   ROMO1   ERCC1   PLA2R1     MIR21   MME     ATM     SERPINE1        WNT16   ATR     TERT    TP53    WRN     CTC1
GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE        http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE   CDKN1A  MAPK14  PLA2R1  SIRT1      ARNTL   WNT16   TP53    MAPKAPK5
GOBP_REGULATION_OF_CELLULAR_SENESCENCE  http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_REGULATION_OF_CELLULAR_SENESCENCE     AKT3    MIR543  CDK6    ZMPSTE24  PLK2     NEK6    ZNF277  CGAS    VASH1   SIRT1   MORC3   ABL1    RSL1D1  FBXO5   HLA-G   ING2    KIR2DL4 ARG2    ARNTL   MIR10A  MIR146A MIR17   MIR188  MIR217  MIR22   MIR34A     MAP3K3  MIF     YBX1    PAWR    ABI3    FZR1    SIRT6   SLC30A10        B2M     PTEN    MIR20B  RBL1    BCL6    BMPR1A  NEK4    TBX2    MIR590  TERC    TERF2   TERT       TP53    TWIST1  WNT1    ZKSCAN3 HMGA2   YPEL3   PNPT1   EEF1E1  NUAK1
GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE    AKT3    MIR543  CDK6       PLK2    SIRT1   ABL1    FBXO5   MIR17   MAP3K3  MIF     YBX1    FZR1    SIRT6   SLC30A10        PTEN    RBL1    BCL6    TBX2    MIR590  TERC    TERF2   TERT    TWIST1     WNT1    ZKSCAN3 HMGA2
GOBP_POSITIVE_REGULATION_OF_CELLULAR_SENESCENCE http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_CELLULAR_SENESCENCE    CGAS    SIRT1   MORC3      HLA-G   KIR2DL4 ARG2    MIR10A  MIR146A MIR188  MIR217  MIR22   MIR34A  PAWR    ABI3    B2M     MIR20B  TP53    YPEL3   EEF1E1
```
