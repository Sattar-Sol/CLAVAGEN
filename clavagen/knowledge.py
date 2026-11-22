"""Small built-in gene-disease knowledge base.

This is a convenience starter set only. It is NOT a substitute for OMIM,
ClinGen gene-disease validity curation, GeneReviews or ClinVar. Users should
supply curated, up-to-date information via an annotation JSON file, which
always overrides the values here. Verify MIM numbers before reporting.
"""

# ACMG SF v3.3 (Genet Med 2025;27:101454): v3.2 list (81 genes) + ABCD1, CYP27A1, PLN.
ACMG_SF_VERSION = "v3.3"
ACMG_SF_GENES = set("""
ABCD1 ACTA2 ACTC1 ACVRL1 APC APOB ATP7B BAG3 BMPR1A BRCA1 BRCA2 BTD CACNA1S CALM1 CALM2 CALM3
CASQ2 COL3A1 CYP27A1 DES DSC2 DSG2 DSP ENG FBN1 FLNC GAA GLA HFE HNF1A KCNH2 KCNQ1 LDLR LMNA
MAX MEN1 MLH1 MSH2 MSH6 MUTYH MYBPC3 MYH11 MYH7 MYL2 MYL3 NF2 OTC PALB2 PCSK9 PKP2 PLN PMS2
PRKAG2 PTEN RB1 RBM20 RET RPE65 RYR1 RYR2 SCN5A SDHAF2 SDHB SDHC SDHD SMAD3 SMAD4 STK11
TGFBR1 TGFBR2 TMEM127 TMEM43 TNNC1 TNNI3 TNNT2 TP53 TPM1 TRDN TSC1 TSC2 TTN TTR VHL WT1
""".split())

GENES = {
    "HBB": {
        "gene": "HBB", "mim": "141900", "inheritance": "AR", "lof_mechanism": True,
        "diseases": [
            {"name": "Sickle cell disease", "mim": "603903", "inheritance": "AR",
             "summary": "Haemoglobinopathy caused by HbS (p.Glu7Val) in homozygous or compound-heterozygous state; "
                        "carriers (sickle cell trait) are usually asymptomatic.",
             "features": ["haemolytic anaemia", "vaso-occlusive pain crises", "acute chest syndrome",
                          "stroke", "splenic sequestration", "increased infection risk"]},
            {"name": "Beta-thalassaemia", "mim": "613985", "inheritance": "AR",
             "summary": "Reduced (beta+) or absent (beta0) beta-globin synthesis; severity ranges from "
                        "carrier (thalassaemia minor) to transfusion-dependent thalassaemia major.",
             "features": ["microcytic hypochromic anaemia", "ineffective erythropoiesis",
                          "hepatosplenomegaly", "iron overload", "bone deformities"]},
        ],
        "phenotype_keywords": ["anemia", "anaemia", "sickle", "thalass", "hemoly", "haemoly", "microcyt",
                               "vaso-occlusive", "splenomegaly", "hbs", "hemoglobin", "haemoglobin"],
        "management": ["Haemoglobin electrophoresis/HPLC to confirm phenotype",
                       "Partner/reproductive carrier testing for AR condition",
                       "Haematology referral for affected individuals"],
    },
    "CFTR": {
        "gene": "CFTR", "mim": "602421", "inheritance": "AR", "lof_mechanism": True,
        "diseases": [{"name": "Cystic fibrosis", "mim": "219700", "inheritance": "AR",
                      "summary": "Multisystem disorder caused by biallelic CFTR dysfunction affecting "
                                 "lungs, pancreas, intestine and reproductive tract.",
                      "features": ["chronic sinopulmonary infection", "bronchiectasis", "pancreatic insufficiency",
                                   "meconium ileus", "elevated sweat chloride", "CBAVD / male infertility"]}],
        "phenotype_keywords": ["cystic", "sweat chloride", "pancrea", "bronchiect", "meconium", "cbavd"],
        "management": ["Sweat chloride testing", "Variant-specific CFTR modulator eligibility review",
                       "CF centre referral"],
    },
    "BRCA1": {
        "gene": "BRCA1", "mim": "113705", "inheritance": "AD", "lof_mechanism": True,
        "diseases": [{"name": "Hereditary breast and ovarian cancer (BRCA1)", "mim": "604370", "inheritance": "AD",
                      "summary": "Cancer predisposition with markedly increased lifetime risk of breast and "
                                 "ovarian cancer; also prostate and pancreatic cancer.",
                      "features": ["early-onset breast cancer", "triple-negative breast cancer", "ovarian cancer",
                                   "pancreatic cancer", "prostate cancer"]}],
        "phenotype_keywords": ["breast", "ovarian", "pancrea", "prostate"],
        "management": ["Refer to cancer genetics", "Enhanced breast surveillance / risk-reducing options",
                       "Cascade testing of relatives", "PARP-inhibitor eligibility in relevant cancers"],
    },
    "BRCA2": {
        "gene": "BRCA2", "mim": "600185", "inheritance": "AD", "lof_mechanism": True,
        "diseases": [{"name": "Hereditary breast and ovarian cancer (BRCA2)", "mim": "612555", "inheritance": "AD",
                      "summary": "Cancer predisposition (breast, ovarian, male breast, prostate, pancreatic). "
                                 "Biallelic variants cause Fanconi anaemia complementation group D1.",
                      "features": ["breast cancer", "male breast cancer", "ovarian cancer", "prostate cancer",
                                   "pancreatic cancer"]}],
        "phenotype_keywords": ["breast", "ovarian", "pancrea", "prostate", "fanconi"],
        "management": ["Refer to cancer genetics", "Cascade testing of relatives"],
    },
    "TP53": {
        "gene": "TP53", "mim": "191170", "inheritance": "AD", "lof_mechanism": True,
        "diseases": [{"name": "Li-Fraumeni syndrome", "mim": "151623", "inheritance": "AD",
                      "summary": "High-penetrance multi-cancer predisposition syndrome.",
                      "features": ["soft-tissue sarcoma", "osteosarcoma", "premenopausal breast cancer",
                                   "brain tumours", "adrenocortical carcinoma", "leukaemia"]}],
        "phenotype_keywords": ["sarcoma", "adrenocortical", "brain tumo", "breast", "choroid plexus", "leuk"],
        "management": ["Toronto-protocol surveillance (incl. whole-body MRI)", "Avoid ionising radiation where possible",
                       "Cascade testing; consider somatic/CHIP origin of low-VAF variants"],
    },
    "MLH1": {
        "gene": "MLH1", "mim": "120436", "inheritance": "AD", "lof_mechanism": True,
        "diseases": [{"name": "Lynch syndrome (MLH1)", "mim": "609310", "inheritance": "AD",
                      "summary": "Mismatch-repair deficiency predisposing to colorectal, endometrial and other cancers.",
                      "features": ["colorectal cancer", "endometrial cancer", "ovarian cancer", "gastric cancer",
                                   "urothelial cancer", "MSI-high tumours"]}],
        "phenotype_keywords": ["colorectal", "colon", "endometri", "msi", "mismatch", "lynch"],
        "management": ["Colonoscopy every 1-2 years", "Tumour MMR IHC / MSI correlation", "Cascade testing"],
    },
    "FBN1": {
        "gene": "FBN1", "mim": "134797", "inheritance": "AD", "lof_mechanism": True,
        "diseases": [{"name": "Marfan syndrome", "mim": "154700", "inheritance": "AD",
                      "summary": "Connective-tissue disorder with cardiovascular, ocular and skeletal involvement.",
                      "features": ["aortic root dilatation/dissection", "ectopia lentis", "tall stature",
                                   "arachnodactyly", "pectus deformity", "mitral valve prolapse"]}],
        "phenotype_keywords": ["aort", "ectopia", "lens", "marfan", "arachnodact", "tall"],
        "management": ["Echocardiography / aortic imaging", "Ophthalmology review", "Cascade testing"],
    },
    "LDLR": {
        "gene": "LDLR", "mim": "606945", "inheritance": "AD (semi-dominant)", "lof_mechanism": True,
        "diseases": [{"name": "Familial hypercholesterolaemia", "mim": "143890", "inheritance": "AD",
                      "summary": "Elevated LDL cholesterol from birth with premature atherosclerotic disease.",
                      "features": ["elevated LDL-C", "tendon xanthomas", "corneal arcus", "premature coronary disease"]}],
        "phenotype_keywords": ["ldl", "cholesterol", "xanthom", "hypercholest", "coronary"],
        "management": ["Lipid-lowering therapy", "Cascade testing (including children)"],
    },
    "PAH": {
        "gene": "PAH", "mim": "612349", "inheritance": "AR", "lof_mechanism": True,
        "diseases": [{"name": "Phenylketonuria / hyperphenylalaninaemia", "mim": "261600", "inheritance": "AR",
                      "summary": "Phenylalanine hydroxylase deficiency; detected on newborn screening.",
                      "features": ["elevated phenylalanine", "intellectual disability if untreated",
                                   "seizures", "eczema", "hypopigmentation"]}],
        "phenotype_keywords": ["phenylalanine", "pku", "hyperphenyl"],
        "management": ["Metabolic referral", "Dietary management", "BH4-responsiveness assessment"],
    },
    "DMD": {
        "gene": "DMD", "mim": "300377", "inheritance": "XLR", "lof_mechanism": True,
        "diseases": [{"name": "Duchenne / Becker muscular dystrophy", "mim": "310200; 300376", "inheritance": "XLR",
                      "summary": "Dystrophinopathy; out-of-frame changes usually cause Duchenne, in-frame often Becker.",
                      "features": ["progressive proximal weakness", "elevated CK", "calf pseudohypertrophy",
                                   "dilated cardiomyopathy"]}],
        "phenotype_keywords": ["muscular dystrophy", "ck", "creatine kinase", "weakness", "cardiomyopathy"],
        "management": ["Neuromuscular referral", "Cardiac surveillance", "Exon-skipping eligibility review",
                       "Carrier testing of mother/female relatives"],
    },
    "SOD1": {
        "gene": "SOD1", "mim": "147450", "inheritance": "AD", "lof_mechanism": False,
        "diseases": [
            {"name": "Amyotrophic lateral sclerosis 1 (ALS1)", "mim": "105400", "inheritance": "AD (rarely AR)",
             "summary": "Adult-onset motor neuron disease; SOD1 variants act mainly through a toxic gain of function "
                        "(protein misfolding/aggregation). Most are dominant; a few (e.g. p.Asp91Ala) are recessive.",
             "features": ["progressive muscle weakness", "upper and lower motor neuron signs", "fasciculations",
                          "bulbar dysfunction", "respiratory failure", "family history of ALS"]},
            {"name": "Spastic tetraplegia and axial hypotonia, progressive (STAHP)", "mim": "618598",
             "inheritance": "AR",
             "summary": "Rare childhood disorder caused by biallelic SOD1 loss of function.",
             "features": ["infantile hypotonia", "progressive spastic tetraplegia"]}],
        "phenotype_keywords": ["als", "amyotrophic", "motor neuron", "fasciculation", "bulbar", "weakness",
                               "spastic", "denervation"],
        "management": ["Neurology / ALS clinic referral", "Assess eligibility for SOD1-targeted therapy (tofersen)",
                       "Genetic counselling; predictive testing of relatives only with counselling (reduced penetrance)"],
    },
    "MYH7": {
        "gene": "MYH7", "mim": "160760", "inheritance": "AD", "lof_mechanism": False,
        "diseases": [{"name": "Hypertrophic cardiomyopathy (MYH7)", "mim": "192600", "inheritance": "AD",
                      "summary": "Sarcomeric cardiomyopathy predominantly caused by missense variants.",
                      "features": ["left ventricular hypertrophy", "arrhythmia", "sudden cardiac death", "heart failure"]}],
        "phenotype_keywords": ["hypertroph", "cardiomyopathy", "lvh", "arrhythm", "sudden"],
        "management": ["Cardiology referral", "Echocardiography/ECG screening of relatives"],
        "missense_constrained": True,
    },
}


def gene_info(symbol):
    if not symbol:
        return {}
    import copy
    return copy.deepcopy(GENES.get(symbol.upper(), {"gene": symbol.upper()}))
