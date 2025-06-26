
AA_COLOUR_SCHEME = {
    'P': 'deeppink',
    'M': 'orange',
    'A': 'orange',
    'V': 'orange',
    'I': 'orange',
    'L': 'orange',
    'F': 'orange',
    'Y': 'orange',
    'W': 'orange',
    'H': 'seagreen',
    'R': 'seagreen',
    'K': 'seagreen',
    'D': 'firebrick',
    'E': 'firebrick',
    'N': 'dodgerblue',
    'Q': 'dodgerblue',
    'S': 'dodgerblue',
    'T': 'dodgerblue',
    'G': 'dodgerblue',
    'C': 'dodgerblue',
    'X': 'black',
}
START_CODONS = [
    'AUG', 'CUG', 'GUG', 'AUC', 'ACG',
]
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AMINO_ACID_GROUPS = {
    'basic': ['K', 'R', 'H'],
    'special-case': ['P'],
    'hydrophobic': ['M', 'A','V','I','L','F','Y','W'],
    'acidic': ['D', 'E'],
    'polar': ['N','Q','S','T','G','C'],
}
BGD_COLOURS = ['wheat','#157F91']

TRAIN_FEATURES = {
    'canonical': [
        'ubiCounts', 'acetylCounts',
        'A_frac', 'C_frac', 'G_frac', 'U_frac',
        'protLength', 'relativePosition', 'start_dist',
        'proteinHydrophobicity', 'peptideHydrophobicity',
        'C_term_acidic', 'C_term_basic', 'C_term_hydrophobic', 'C_term_polar', 'C_term_special-case',
        'C_term_neg_1_acidic', 'C_term_neg_1_basic', 'C_term_neg_1_end',
        'C_term_neg_1_hydrophobic', 'C_term_neg_1_polar', 'C_term_neg_1_special-case',
    ],
    'cryptic': [
        'stratum',
    
        'protLength', 'relativePosition', 'proteinHydrophobicity', 'peptideHydrophobicity',
    
        'C_term_basic', 'C_term_special-case', 'C_term_hydrophobic', 'C_term_acidic', 'C_term_polar',
        
        'C_term_neg_1_basic', 'C_term_neg_1_special-case', 'C_term_neg_1_hydrophobic', 'C_term_neg_1_acidic', 'C_term_neg_1_polar', 'C_term_neg_1_end',

        'A_frac', 'C_frac', 'G_frac', 'U_frac',

        'codingStrand', 'mismatches', 'deNovoAssembly',

        'stopDistances', 'start_dist', 'kozakScore',
    ] + [f'{start_codon}_upstream' for start_codon in START_CODONS],
    'spliced': (
        [
            'ubiCounts', 'acetylCounts',
            # 'UTR_length', 'UTR_A_frac', 'UTR_C_frac', 'UTR_G_frac', 'UTR_U_frac',
            'nCanonicalPeptides', 'sr2_can_dist', 'sr1_can_dist',
            'interveningSeqLengths',
            'isForward',
            'protLength',
            # 'Cytoplasmic translation', 'RNA processing', 'Nucleic acid metabolic process',
        ] + [
            f'p1_{amino_acid}' for amino_acid in AMINO_ACIDS
        ]  
        + [
            f'p1_prime_{amino_acid}' for amino_acid in AMINO_ACIDS
        ] 
        # + [
        #     f'p_neg_1_{amino_acid}' for amino_acid in AMINO_ACIDS
        # ] + [
        #     f'p_neg_1_prime_{amino_acid}' for amino_acid in AMINO_ACIDS
        # ]
    )

}


COLOUR_DICT = {
    'canonical': '#EC9A56',
    'spliced': '#9BBFE5',
    'cryptic': '#BA69BE',
}
COMMON_FEATURES = [
    'proteinID', 'protSeq', 'rnaSeq', 'iupred3_preds', 'proteinHydrophobicity'
]
CRYPTIC_STRATA = [
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
]


STOP_CODONS = [
    'UAA', 'UAG', 'UGA'
]
N_CV_GROUPS = 10

TRANSCRIPT_FEATURES = {
    'K562': ['tr_TPM_K562_bulk', 'tr_TPM_K562_free', 'tr_TPM_K562_S80', 'tr_TPM_K562_poly'],
    'B721.221': ['tr_TPM_721'],
}
PROTEOMICS_FEATURES = {
    'K562': ['proteomics_K562'],
    'B721.221': ['proteomics_B721'],
}
ID_COLUMNS = ['peptide', 'proteinID']
NUCLEOTIDE_COLOUR_SCHEME = {
    'A': '#FFAA33',
    'C': 'darkgreen',
    'G': '#1f618d',
    'U': '#FA2A55',
}

STRATUM_SPECIFIC_FEATURES = {
    'canonical': [
        'geneID',
        'ubiCounts', 'acetylCounts',
        # 'UTR_A_frac', 'UTR_C_frac', 'UTR_G_frac', 'UTR_U_frac',
    ],
    'spliced': ['geneID', 'nCanonicalPeptides', 'canonicalPeptides', 'ubiCounts', 'acetylCounts', 'Cytoplasmic translation', 'RNA processing', 'Nucleic acid metabolic process'],
    'fiveUTR': ['geneID',],
    'threeUTR': ['geneID',],
    'CDS_frameshift': ['geneID', 'codingStrand', 'mismatches'],
    'lncRNA': [],
    'intronic': ['geneID'],
    'intergenic': ['deNovoAssembly'],
}


STRATUM_COLOUR_SCHEME = {
    'fiveUTR': 'orange',
    'threeUTR': 'yellow',
    'CDS_frameshift': 'purple',
    'lncRNA': 'darkgrey',
    'intronic': 'forestgreen',
    'intergenic': 'navy',
}

SPLICE_SPECIFIC_FEATURES = [ 'sr1', 'interveningSeqLengths', 'sr1_Index', 'sr2_Index', 'isForward']
SPLICED_FEATURES = (
    [
        'ubiCounts', 'acetylCounts',
        'nCanonicalPeptides', 'sr2_can_dist', 'sr1_can_dist',
        'interveningSeqLengths',
        'isForward',
        'protLength', #'SLIDER_score',
        # 'sr1_localDisorder', 'sr2_localDisorder',
    ] + [
        f'p1_{amino_acid}' for amino_acid in AMINO_ACIDS
    ]  
    + [
        f'p1_prime_{amino_acid}' for amino_acid in AMINO_ACIDS
    ] 
    # + [
    #     f'p_neg_1_{amino_acid}' for amino_acid in AMINO_ACIDS
    # ]  + [
    #     f'p_neg_1_prime_{amino_acid}' for amino_acid in AMINO_ACIDS
    # ]
)
