import pandas as pd
import numpy as np
import xlsxwriter

### STRUCTURAL VARIANTS IN OPERON ###

# Import files
df_ch1 = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/rnaseq_atcc/ATCC23344_v2_CH1_OperonSEQer_voteThreshold1_operonList.csv', header=None)
df_ch2 = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/rnaseq_atcc/ATCC23344_v2_CH2_OperonSEQer_voteThreshold1_operonList.csv', header=None)
gff_ch1 = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/OperonSEQer/23344_CH1_gene.gff', header=None, delimiter='\t')
gff_ch2 = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/OperonSEQer/23344_CH2_gene.gff', header=None, delimiter='\t')
gff_all = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/OperonSEQer/23344_CH1_CH2.gff', header=None, delimiter='\t')
df_rnaseq = pd.read_csv('/media/Storage/shared/Glanders/RNA_SEQ/test_table_new_1.csv')

# Clean-up operon files
df_ch1 = df_ch1.replace('-','_', regex=True)
df_ch2 = df_ch2.replace('-','_', regex=True)
df_ch1 = df_ch1.rename(columns={0: 'GENES'})
df_ch2 = df_ch2.rename(columns={0: 'GENES'})

# Clean-up GFF file
gff_all = gff_all[gff_all[2].str.contains('CDS', na=False)]
gff_all['WP_ID'] = gff_all[8].str.split('ID=cds-WP_', expand = True)[1].str.split(';', expand = True)[0].str.split('-', expand = True)[0]
gff_all['BMA_ID'] = gff_all[8].str.split('Parent=gene-', expand = True)[1].str.split(';', expand = True)[0]
gff_all['GENE'] = gff_all[8].str.split(';gene=', expand = True)[1].str.split(';', expand = True)[0]

# Clean-up RNAseq file
df_rnaseq['BMA_ID'] = np.where(df_rnaseq['Locus_tag.gene'].str.contains('BMA_RS'), df_rnaseq['Locus_tag.gene'], None)
df_rnaseq_BMA =  df_rnaseq[df_rnaseq['BMA_ID'].str.contains('BMA_RS', na=False)].reset_index(drop = True)
df_rnaseq_none = df_rnaseq[~df_rnaseq['BMA_ID'].str.contains('BMA_RS', na=False)].reset_index(drop = True)

df_rnaseq_none['BMA_ID'] = np.where(df_rnaseq_none['target_id'].str.contains('cds_WP'), df_rnaseq_none['target_id'], None)
df_rnaseq_WP = df_rnaseq_none[df_rnaseq_none['BMA_ID'].str.contains('cds_WP', na=False)].reset_index(drop = True)
df_rnaseq_none = df_rnaseq_none[~df_rnaseq_none['BMA_ID'].str.contains('cds_WP', na=False)].reset_index(drop = True)

df_rnaseq_WP['BMA_ID'] = df_rnaseq_WP.BMA_ID.str.split('_', expand = True)[4]
for i in range(len(df_rnaseq_WP)):
    for j in range(len(gff_all)):
        if df_rnaseq_WP.iloc[i]['BMA_ID'] == gff_all.iloc[j]['WP_ID']:
            df_rnaseq_WP.at[i, 'BMA_ID'] = gff_all.iloc[j]['BMA_ID']
            break

for i in range(len(df_rnaseq_none)):
    for j in range(len(gff_all)):
        if df_rnaseq_none.iloc[i]['Locus_tag.gene'] == gff_all.iloc[j]['GENE']:
            df_rnaseq_none.at[i, 'BMA_ID'] = gff_all.iloc[j]['BMA_ID']
            break

df_rnaseq_final = pd.concat([df_rnaseq_BMA, df_rnaseq_WP], axis = 0, ignore_index=True)
df_rnaseq_final = pd.concat([df_rnaseq_final, df_rnaseq_none], axis = 0, ignore_index=True)
df_rnaseq_final = df_rnaseq_final.sort_values(by=['Unnamed: 0'], ignore_index=True)
df_rnaseq_final['Locus_tag.gene'] = df_rnaseq_final['BMA_ID']
df_rnaseq_final = df_rnaseq_final.drop(['BMA_ID'], axis=1)
df_rnaseq_filt = df_rnaseq_final[(df_rnaseq_final['qval'] < 1E-3) & ((df_rnaseq_final['log2FC'] < -1) | (df_rnaseq_final['log2FC'] > 1))].reset_index(drop = True)

# Generate operon information using GFF files for chromosome 1
df_ch1['START_GENE'] = df_ch1['GENES'].str.split(';').str[0]
df_ch1['END_GENE'] = df_ch1['GENES'].str.split(';').str[-1]
df_ch1['CHROM'] = None
df_ch1['START'] = None
df_ch1['END'] = None
df_ch1['STRAND'] = None
for i in range(len(df_ch1)):
    df_ch1.loc[i,'CHROM'] = gff_ch1[gff_ch1[8].str.contains(df_ch1.iloc[i]['START_GENE'])][0].to_string(index=False)
    df_ch1.loc[i,'START'] = gff_ch1[gff_ch1[8].str.contains(df_ch1.iloc[i]['START_GENE'])][3].to_string(index=False)
    df_ch1.loc[i,'END'] = gff_ch1[gff_ch1[8].str.contains(df_ch1.iloc[i]['END_GENE'])][4].to_string(index=False)
    df_ch1.loc[i,'STRAND'] = gff_ch1[gff_ch1[8].str.contains(df_ch1.iloc[i]['START_GENE'])][6].to_string(index=False)

df_ch1['STRAND'] = df_ch1['STRAND'].str.replace('+','1', regex=False)
df_ch1['STRAND'] = df_ch1['STRAND'].str.replace('-','-1', regex=False)
df_ch1 = df_ch1.drop(['START_GENE', 'END_GENE'], axis=1)
df_ch1['OPERON_ID'] = None
for i in range(len(df_ch1)):
    df_ch1.loc[i,'OPERON_ID'] = "OPERON_{0:04}".format(i+1)
    
# Generate operon information using GFF files for chromosome 2
df_ch2['START_GENE'] = df_ch2['GENES'].str.split(';').str[0]
df_ch2['END_GENE'] = df_ch2['GENES'].str.split(';').str[-1]
df_ch2['CHROM'] = None
df_ch2['START'] = None
df_ch2['END'] = None
df_ch2['STRAND'] = None
for i in range(len(df_ch2)):
    df_ch2.loc[i,'CHROM'] = gff_ch2[gff_ch2[8].str.contains(df_ch2.iloc[i]['START_GENE'])][0].to_string(index=False)
    df_ch2.loc[i,'START'] = gff_ch2[gff_ch2[8].str.contains(df_ch2.iloc[i]['START_GENE'])][3].to_string(index=False)
    df_ch2.loc[i,'END'] = gff_ch2[gff_ch2[8].str.contains(df_ch2.iloc[i]['END_GENE'])][4].to_string(index=False)
    df_ch2.loc[i,'STRAND'] = gff_ch2[gff_ch2[8].str.contains(df_ch2.iloc[i]['START_GENE'])][6].to_string(index=False)

df_ch2['STRAND'] = df_ch2['STRAND'].str.replace('+','1', regex=False)
df_ch2['STRAND'] = df_ch2['STRAND'].str.replace('-','-1', regex=False)
df_ch2 = df_ch2.drop(['START_GENE', 'END_GENE'], axis=1)
df_ch2['OPERON_ID'] = None
for i in range(len(df_ch2)):
    df_ch2.loc[i,'OPERON_ID'] = "OPERON_{0:04}".format(i+1+len(df_ch1))

# Find operons disrupted by structural variants
disrupted_operon_df_list = list()
disrupted_operon_df_all = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', 'SV_TYPE'])
summary_operon = pd.DataFrame(columns=['SAMPLE', '# SVS', 'DISRUPTED OPERONS', 'DISRUPTED GENES'])
tmp1 = df_ch1.assign(GENE=df_ch1.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
tmp2 = df_ch2.assign(GENE=df_ch2.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
total_genes = len(tmp1) + len(tmp2)
sample_list = list(pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/mumandco_v2/summary.all.sv.csv').index)
for sample in sample_list:
    try:
        sv_list = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/mumandco_v2/' + sample + '_output/' + sample + '.SVs_all.tsv', sep="\t")
        sv_list_ch1 = sv_list[sv_list['ref_chr'] == 'NC_006348.1']
        sv_list_ch2 = sv_list[sv_list['ref_chr'] == 'NC_006349.2']
        sv_list_ch1_del_contr = sv_list_ch1[(sv_list_ch1['SV_type'] == 'contraction') | (sv_list_ch1['SV_type'] == 'deletion_novel') | (sv_list_ch1['SV_type'] == 'deletion_mobile') | (sv_list_ch1['SV_type'] == 'deletion_manual')].reset_index(drop = True)
        sv_list_ch1_other = sv_list_ch1[~((sv_list_ch1['SV_type'] == 'contraction') | (sv_list_ch1['SV_type'] == 'deletion_novel') | (sv_list_ch1['SV_type'] == 'deletion_mobile') | (sv_list_ch1['SV_type'] == 'deletion_manual'))].reset_index(drop = True)
        sv_list_ch2_del_contr = sv_list_ch2[(sv_list_ch2['SV_type'] == 'contraction') | (sv_list_ch2['SV_type'] == 'deletion_novel') | (sv_list_ch2['SV_type'] == 'deletion_mobile') | (sv_list_ch2['SV_type'] == 'deletion_manual')].reset_index(drop = True)
        sv_list_ch2_other = sv_list_ch2[~((sv_list_ch2['SV_type'] == 'contraction') | (sv_list_ch2['SV_type'] == 'deletion_novel') | (sv_list_ch2['SV_type'] == 'deletion_mobile') | (sv_list_ch2['SV_type'] == 'deletion_manual'))].reset_index(drop = True)

        disrupted_operon_df = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', 'SV_TYPE'])

        for i in range(len(df_ch1)):
            for j in range(len(sv_list_ch1_other)):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1.iloc[i]['OPERON_ID']],
                    'START': [df_ch1.iloc[i]['START']],
                    'END': [df_ch1.iloc[i]['END']],
                    'STRAND': [df_ch1.iloc[i]['STRAND']],
                    'GENES': [df_ch1.iloc[i]['GENES']],
                    'SV_TYPE': [sv_list_ch1_other.iloc[j]['SV_type']]})
                if (int(df_ch1.iloc[i]['START']) < sv_list_ch1_other.iloc[j]['ref_start']) & (int(df_ch1.iloc[i]['END']) > sv_list_ch1_other.iloc[j]['ref_start']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch1.iloc[i]['START']) < sv_list_ch1_other.iloc[j]['ref_stop']) & (int(df_ch1.iloc[i]['END']) > sv_list_ch1_other.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch1.iloc[i]['START']) > sv_list_ch1_other.iloc[j]['ref_start']) & (int(df_ch1.iloc[i]['START']) < sv_list_ch1_other.iloc[j]['ref_start'] + 5000):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch1.iloc[i]['END']) < sv_list_ch1_other.iloc[j]['ref_stop']) & (int(df_ch1.iloc[i]['START']) > sv_list_ch1_other.iloc[j]['ref_stop'] - 5000):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)

            for j in range(len(sv_list_ch1_del_contr)):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1.iloc[i]['OPERON_ID']],
                    'START': [df_ch1.iloc[i]['START']],
                    'END': [df_ch1.iloc[i]['END']],
                    'STRAND': [df_ch1.iloc[i]['STRAND']],
                    'GENES': [df_ch1.iloc[i]['GENES']],
                    'SV_TYPE': [sv_list_ch1_del_contr.iloc[j]['SV_type']]})
                if (int(df_ch1.iloc[i]['START']) < sv_list_ch1_del_contr.iloc[j]['ref_start']) & (int(df_ch1.iloc[i]['END']) > sv_list_ch1_del_contr.iloc[j]['ref_start']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch1.iloc[i]['START']) > sv_list_ch1_del_contr.iloc[j]['ref_start']) & (int(df_ch1.iloc[i]['END']) < sv_list_ch1_del_contr.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch1.iloc[i]['START']) < sv_list_ch1_del_contr.iloc[j]['ref_stop']) & (int(df_ch1.iloc[i]['END']) > sv_list_ch1_del_contr.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)

        for i in range(len(df_ch2)):
            for j in range(len(sv_list_ch2_other)):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2.iloc[i]['OPERON_ID']],
                    'START': [df_ch2.iloc[i]['START']],
                    'END': [df_ch2.iloc[i]['END']],
                    'STRAND': [df_ch2.iloc[i]['STRAND']],
                    'GENES': [df_ch2.iloc[i]['GENES']],
                    'SV_TYPE': [sv_list_ch2_other.iloc[j]['SV_type']]})
                if (int(df_ch2.iloc[i]['START']) < sv_list_ch2_other.iloc[j]['ref_start']) & (int(df_ch2.iloc[i]['END']) > sv_list_ch2_other.iloc[j]['ref_start']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch2.iloc[i]['START']) < sv_list_ch2_other.iloc[j]['ref_stop']) & (int(df_ch2.iloc[i]['END']) > sv_list_ch2_other.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch2.iloc[i]['START']) > sv_list_ch2_other.iloc[j]['ref_start']) & (int(df_ch2.iloc[i]['START']) < sv_list_ch2_other.iloc[j]['ref_start'] + 5000):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch2.iloc[i]['END']) < sv_list_ch2_other.iloc[j]['ref_stop']) & (int(df_ch2.iloc[i]['START']) > sv_list_ch2_other.iloc[j]['ref_stop'] - 5000):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
            for j in range(len(sv_list_ch2_del_contr)):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2.iloc[i]['OPERON_ID']],
                    'START': [df_ch2.iloc[i]['START']],
                    'END': [df_ch2.iloc[i]['END']],
                    'STRAND': [df_ch2.iloc[i]['STRAND']],
                    'GENES': [df_ch2.iloc[i]['GENES']],
                    'SV_TYPE': [sv_list_ch2_del_contr.iloc[j]['SV_type']]})
                if (int(df_ch2.iloc[i]['START']) < sv_list_ch2_del_contr.iloc[j]['ref_start']) & (int(df_ch2.iloc[i]['END']) > sv_list_ch2_del_contr.iloc[j]['ref_start']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch2.iloc[i]['START']) > sv_list_ch2_del_contr.iloc[j]['ref_start']) & (int(df_ch2.iloc[i]['END']) < sv_list_ch2_del_contr.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)
                elif (int(df_ch2.iloc[i]['START']) < sv_list_ch2_del_contr.iloc[j]['ref_stop']) & (int(df_ch2.iloc[i]['END'])> sv_list_ch2_del_contr.iloc[j]['ref_stop']):
                    disrupted_operon_df = pd.concat([disrupted_operon_df, new_entry], ignore_index=True)

        disrupted_operon_df['SV_TYPE'] = disrupted_operon_df.groupby(['CHROM', 'START'])['SV_TYPE'].transform(lambda x: '; '.join(x))
        disrupted_operon_df = disrupted_operon_df.drop_duplicates().reset_index(drop = True)
    except:
        disrupted_operon_df = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', 'SV_TYPE'])
        sv_list = []
    

    disrupted_operon_df_list.append(disrupted_operon_df)
    disrupted_operon_df_all = pd.concat([disrupted_operon_df_all, disrupted_operon_df], ignore_index=True)
                                                             
    n_operon = len(df_ch1) + len(df_ch2)
    n_sv = len(sv_list)
    n_broken_operon = len(disrupted_operon_df)
    p_broken_operon = round(n_broken_operon * 100 / n_operon, 1)
    broken_operon = str(n_broken_operon) + " ("+ str(p_broken_operon) + "%)"
    
    disrupted_operon_df = disrupted_operon_df.assign(GENE=disrupted_operon_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
    n_genes = len(disrupted_operon_df)
    p_genes = round(n_genes * 100 / total_genes, 1)
    broken_genes = str(n_genes) + " ("+ str(p_genes) + "%)" 
    
    new_summary = pd.DataFrame.from_dict({
                        'SAMPLE': [sample],
                        '# SVS': [n_sv],
                        'DISRUPTED OPERONS': [broken_operon],
                        'DISRUPTED GENES': [broken_genes]
    })
    summary_operon = pd.concat([summary_operon, new_summary], ignore_index=True)


# Filter operons to get Zagreb specific results
zagreb_operon_df = disrupted_operon_df_all[disrupted_operon_df_all['SAMPLE'].str.contains('Zagreb', na=False)].reset_index(drop = True)
zagreb_operon_df = zagreb_operon_df.assign(GENE=zagreb_operon_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
zagreb_operon_df = pd.merge(zagreb_operon_df, df_rnaseq_final[['Locus_tag.gene', 'log2FC', 'qval']], left_on='GENE', right_on='Locus_tag.gene', how='left').reset_index(drop = True)
zagreb_operon_df = zagreb_operon_df.drop(['GENE', 'Locus_tag.gene'], axis=1)
zagreb_operon_df['log2FC_str'] = round(zagreb_operon_df['log2FC'],3).astype(str)
zagreb_operon_df['LOG2FC'] = zagreb_operon_df.groupby(['OPERON_ID'])['log2FC_str'].transform(lambda x: ';'.join(x))
zagreb_operon_df = zagreb_operon_df[zagreb_operon_df['qval'] < 1E-3]
zagreb_operon_df = zagreb_operon_df[abs(zagreb_operon_df['log2FC']) > 1]
zagreb_operon_df = zagreb_operon_df.drop(['log2FC', 'log2FC_str', 'qval'], axis=1).drop_duplicates().reset_index(drop = True)

# Write results to as xlsx sheets and csv files
writer = pd.ExcelWriter('operon_sv_sheets_V3.xlsx', engine='xlsxwriter')
for i in range(len(disrupted_operon_df_list)):
    disrupted_operon_df_list[i].to_excel(writer, sheet_name=sample_list[i], index=False)
writer.close()

summary_operon.to_csv('operon_sv_summary_V3.csv', index=False)
disrupted_operon_df_all.to_csv('operon_sv_all_V3.csv', index=False)
zagreb_operon_df.to_csv('zagreb_sv_rna_V3.csv', index=False)


### INDELS IN "PROMOTER" REGION OF OPERONS ###

df_ch1['START'] = df_ch1['START'].astype(int)
df_ch2['START'] = df_ch2['START'].astype(int)
df_ch1['END'] = df_ch1['END'].astype(int)
df_ch2['END'] = df_ch2['END'].astype(int)

df_ch1_pos = df_ch1[df_ch1['STRAND'] == "1"].reset_index(drop = True)
df_ch1_neg = df_ch1[df_ch1['STRAND'] == "-1"].reset_index(drop = True)
df_ch2_pos = df_ch2[df_ch2['STRAND'] == "1"].reset_index(drop = True)
df_ch2_neg = df_ch2[df_ch2['STRAND'] == "-1"].reset_index(drop = True)

# Find promoter regions disrupted by indels
disrupted_promoter_df_list = list()
summary_promoter_indel = pd.DataFrame(columns=['SAMPLE', '# INDELS', 'DISRUPTED PROMOTERS', 'DISRUPTED GENES'])
disrupted_promoter_df_all = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', '# INDELS'])

for sample in sample_list:
    variant_list = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/snippy_v2/' + sample + '/snps.csv')
    indel_list = variant_list[variant_list['REF'].str.len() != variant_list['ALT'].str.len()].copy()
    indel_list['STOP_R'] = indel_list['POS'] + indel_list['REF'].str.len()
    indel_list['STOP_A'] = indel_list['POS'] + indel_list['ALT'].str.len()
    indel_list['STOP'] = indel_list[['STOP_R', 'STOP_A']].max(axis=1)
    indel_list = indel_list.drop(['STOP_R', 'STOP_A'], axis=1)
    indel_list_ch1 = indel_list[indel_list['CHROM'] == 'NC_006348.1'].reset_index(drop = True)
    indel_list_ch2 = indel_list[indel_list['CHROM'] == 'NC_006349.2'].reset_index(drop = True)

    disrupted_promoter_df = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', '# INDELS'])
    for j in range(len(indel_list_ch1)):
        for i in range(len(df_ch1_pos)):
            if (((indel_list_ch1.iloc[j]['POS'] > (df_ch1_pos.iloc[i]['START'] - 300)) & (indel_list_ch1.iloc[j]['POS'] < (df_ch1_pos.iloc[i]['START']))) | ((indel_list_ch1.iloc[j]['STOP'] > (df_ch1_pos.iloc[i]['START'] - 300)) & (indel_list_ch1.iloc[j]['STOP'] < (df_ch1_pos.iloc[i]['START'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1_pos.iloc[i]['OPERON_ID']],
                    'START': [df_ch1_pos.iloc[i]['START']],
                    'END': [df_ch1_pos.iloc[i]['END']],
                    'STRAND': [df_ch1_pos.iloc[i]['STRAND']],
                    'GENES': [df_ch1_pos.iloc[i]['GENES']],
                    '# INDELS': [indel_list_ch1.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
        for i in range(len(df_ch1_neg)):
            if (((indel_list_ch1.iloc[j]['POS'] > (df_ch1_neg.iloc[i]['END'])) & (indel_list_ch1.iloc[j]['POS'] < (df_ch1_neg.iloc[i]['END'] + 300))) | ((indel_list_ch1.iloc[j]['STOP'] > (df_ch1_neg.iloc[i]['END'])) & (indel_list_ch1.iloc[j]['STOP'] < (df_ch1_neg.iloc[i]['END'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1_neg.iloc[i]['OPERON_ID']],
                    'START': [df_ch1_neg.iloc[i]['START']],
                    'END': [df_ch1_neg.iloc[i]['END']],
                    'STRAND': [df_ch1_neg.iloc[i]['STRAND']],
                    'GENES': [df_ch1_neg.iloc[i]['GENES']],
                    '# INDELS': [indel_list_ch1.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
    for j in range(len(indel_list_ch2)):
        for i in range(len(df_ch2_pos)):
            if (((indel_list_ch2.iloc[j]['POS'] > (df_ch2_pos.iloc[i]['START'] - 300)) & (indel_list_ch2.iloc[j]['POS'] < (df_ch2_pos.iloc[i]['START']))) | ((indel_list_ch2.iloc[j]['STOP'] > (df_ch2_pos.iloc[i]['START'] - 300)) & (indel_list_ch2.iloc[j]['STOP'] < (df_ch2_pos.iloc[i]['START'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2_pos.iloc[i]['OPERON_ID']],
                    'START': [df_ch2_pos.iloc[i]['START']],
                    'END': [df_ch2_pos.iloc[i]['END']],
                    'STRAND': [df_ch2_pos.iloc[i]['STRAND']],
                    'GENES': [df_ch2_pos.iloc[i]['GENES']],
                    '# INDELS': [indel_list_ch2.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
        for i in range(len(df_ch2_neg)):
            if (((indel_list_ch2.iloc[j]['POS'] > (df_ch2_neg.iloc[i]['END'])) & (indel_list_ch2.iloc[j]['POS'] < (df_ch2_neg.iloc[i]['END'] + 300))) | ((indel_list_ch2.iloc[j]['STOP'] > (df_ch2_neg.iloc[i]['END'])) & (indel_list_ch2.iloc[j]['STOP'] < (df_ch2_neg.iloc[i]['END'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2_neg.iloc[i]['OPERON_ID']],
                    'START': [df_ch2_neg.iloc[i]['START']],
                    'END': [df_ch2_neg.iloc[i]['END']],
                    'STRAND': [df_ch2_neg.iloc[i]['STRAND']],
                    'GENES': [df_ch2_neg.iloc[i]['GENES']],
                    '# INDELS': [indel_list_ch2.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
    disrupted_promoter_df['# INDELS'] = disrupted_promoter_df.groupby(['OPERON_ID'])['# INDELS'].transform('count')
    disrupted_promoter_df = disrupted_promoter_df.drop_duplicates().reset_index(drop = True)
    
    disrupted_promoter_df_list.append(disrupted_promoter_df)
    disrupted_promoter_df_all = pd.concat([disrupted_promoter_df_all, disrupted_promoter_df], ignore_index=True)

    n_operon = len(df_ch1) + len(df_ch2)
    n_indel = len(indel_list)
    n_broken_promoter = len(disrupted_promoter_df)
    p_broken_promoter = round(n_broken_promoter * 100 / n_operon, 1)
    broken_promoter = str(n_broken_promoter) + " ("+ str(p_broken_promoter) + "%)"

    disrupted_promoter_df = disrupted_promoter_df.assign(GENE=disrupted_promoter_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
    n_genes = len(disrupted_promoter_df)
    p_genes = round(n_genes * 100 / total_genes, 1)
    broken_genes = str(n_genes) + " ("+ str(p_genes) + "%)" 
    
    new_summary = pd.DataFrame.from_dict({
                        'SAMPLE': [sample],
                        '# INDELS': [n_indel],
                        'DISRUPTED PROMOTERS': [broken_promoter],
                        'DISRUPTED GENES': [broken_genes]
    })
    summary_promoter_indel = pd.concat([summary_promoter_indel, new_summary], ignore_index=True)

# Filter promoters to get Zagreb specific results
zagreb_promoter_df = disrupted_promoter_df_all[disrupted_promoter_df_all['SAMPLE'].str.contains('Zagreb', na=False)].reset_index(drop = True)
zagreb_promoter_df = zagreb_promoter_df.assign(GENE=zagreb_promoter_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
zagreb_promoter_df = pd.merge(zagreb_promoter_df, df_rnaseq_final[['Locus_tag.gene', 'log2FC', 'qval']], left_on='GENE', right_on='Locus_tag.gene', how='left').reset_index(drop = True)
zagreb_promoter_df = zagreb_promoter_df.drop(['GENE', 'Locus_tag.gene'], axis=1)
zagreb_promoter_df['log2FC_str'] = round(zagreb_promoter_df['log2FC'],3).astype(str)
zagreb_promoter_df['LOG2FC'] = zagreb_promoter_df.groupby(['OPERON_ID'])['log2FC_str'].transform(lambda x: ';'.join(x))
zagreb_promoter_df = zagreb_promoter_df[zagreb_promoter_df['qval'] < 1E-3]
zagreb_promoter_df = zagreb_promoter_df[abs(zagreb_promoter_df['log2FC']) > 1]
zagreb_promoter_df = zagreb_promoter_df.drop(['log2FC', 'log2FC_str', 'qval'], axis=1).drop_duplicates().reset_index(drop = True)

# Write results to as xlsx sheets and csv files
writer = pd.ExcelWriter('promoter_indel_sheets_V3.xlsx', engine='xlsxwriter')
for i in range(len(disrupted_promoter_df_list)):
    disrupted_promoter_df_list[i].to_excel(writer, sheet_name=sample_list[i], index=False)
writer.close()

summary_promoter_indel.to_csv('promoter_indel_summary_V3.csv', index=False)
disrupted_promoter_df_all.to_csv('promoter_indel_all_V3.csv', index=False)
zagreb_promoter_df.to_csv('zagreb_promoter_indel_rna_V3.csv', index=False)


### SNPS IN "PROMOTER" REGION OF OPERONS ###

# Find promoter regions disrupted by snps
disrupted_promoter_snp_df_list = list()
summary_promoter_snp = pd.DataFrame(columns=['SAMPLE', '# SNPS', 'DISRUPTED PROMOTERS', 'DISRUPTED GENES'])
disrupted_promoter_snp_df_all = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', '# SNPS'])

for sample in sample_list:
    print(sample)
    variant_list = pd.read_csv('/media/Data/phil/RedmineIssues/Open/Support27140/snippy_v2/' + sample + '/snps.csv')
    snp_list = variant_list[variant_list['REF'].str.len() == variant_list['ALT'].str.len()].copy()
    snp_list['STOP_R'] = snp_list['POS'] + snp_list['REF'].str.len()
    snp_list['STOP_A'] = snp_list['POS'] + snp_list['ALT'].str.len()
    snp_list['STOP'] = snp_list[['STOP_R', 'STOP_A']].max(axis=1)
    snp_list = snp_list.drop(['STOP_R', 'STOP_A'], axis=1)
    snp_list_ch1 = snp_list[snp_list['CHROM'] == 'NC_006348.1'].reset_index(drop = True)
    snp_list_ch2 = snp_list[snp_list['CHROM'] == 'NC_006349.2'].reset_index(drop = True)

    disrupted_promoter_df = pd.DataFrame(columns=['SAMPLE', 'CHROM', 'OPERON_ID', 'START', 'END', 'STRAND', 'GENES', '# SNPS'])
    for j in range(len(snp_list_ch1)):
        for i in range(len(df_ch1_pos)):
            if (((snp_list_ch1.iloc[j]['POS'] > (df_ch1_pos.iloc[i]['START'] - 300)) & (snp_list_ch1.iloc[j]['POS'] < (df_ch1_pos.iloc[i]['START']))) | ((snp_list_ch1.iloc[j]['STOP'] > (df_ch1_pos.iloc[i]['START'] - 300)) & (snp_list_ch1.iloc[j]['STOP'] < (df_ch1_pos.iloc[i]['START'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1_pos.iloc[i]['OPERON_ID']],
                    'START': [df_ch1_pos.iloc[i]['START']],
                    'END': [df_ch1_pos.iloc[i]['END']],
                    'STRAND': [df_ch1_pos.iloc[i]['STRAND']],
                    'GENES': [df_ch1_pos.iloc[i]['GENES']],
                    '# SNPS': [snp_list_ch1.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
        for i in range(len(df_ch1_neg)):
            if (((snp_list_ch1.iloc[j]['POS'] > (df_ch1_neg.iloc[i]['END'])) & (snp_list_ch1.iloc[j]['POS'] < (df_ch1_neg.iloc[i]['END'] + 300))) | ((snp_list_ch1.iloc[j]['STOP'] > (df_ch1_neg.iloc[i]['END'])) & (snp_list_ch1.iloc[j]['STOP'] < (df_ch1_neg.iloc[i]['END'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006348.1'],
                    'OPERON_ID': [df_ch1_neg.iloc[i]['OPERON_ID']],
                    'START': [df_ch1_neg.iloc[i]['START']],
                    'END': [df_ch1_neg.iloc[i]['END']],
                    'STRAND': [df_ch1_neg.iloc[i]['STRAND']],
                    'GENES': [df_ch1_neg.iloc[i]['GENES']],
                    '# SNPS': [snp_list_ch1.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
    for j in range(len(snp_list_ch2)):
        for i in range(len(df_ch2_pos)):
            if (((snp_list_ch2.iloc[j]['POS'] > (df_ch2_pos.iloc[i]['START'] - 300)) & (snp_list_ch2.iloc[j]['POS'] < (df_ch2_pos.iloc[i]['START']))) | ((snp_list_ch2.iloc[j]['STOP'] > (df_ch2_pos.iloc[i]['START'] - 300)) & (snp_list_ch2.iloc[j]['STOP'] < (df_ch2_pos.iloc[i]['START'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2_pos.iloc[i]['OPERON_ID']],
                    'START': [df_ch2_pos.iloc[i]['START']],
                    'END': [df_ch2_pos.iloc[i]['END']],
                    'STRAND': [df_ch2_pos.iloc[i]['STRAND']],
                    'GENES': [df_ch2_pos.iloc[i]['GENES']],
                    '# SNPS': [snp_list_ch2.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
        for i in range(len(df_ch2_neg)):
            if (((snp_list_ch2.iloc[j]['POS'] > (df_ch2_neg.iloc[i]['END'])) & (snp_list_ch2.iloc[j]['POS'] < (df_ch2_neg.iloc[i]['END'] + 300))) | ((snp_list_ch2.iloc[j]['STOP'] > (df_ch2_neg.iloc[i]['END'])) & (snp_list_ch2.iloc[j]['STOP'] < (df_ch2_neg.iloc[i]['END'])))):
                new_entry = pd.DataFrame.from_dict({
                    'SAMPLE': [sample],
                    'CHROM': ['NC_006349.2'],
                    'OPERON_ID': [df_ch2_neg.iloc[i]['OPERON_ID']],
                    'START': [df_ch2_neg.iloc[i]['START']],
                    'END': [df_ch2_neg.iloc[i]['END']],
                    'STRAND': [df_ch2_neg.iloc[i]['STRAND']],
                    'GENES': [df_ch2_neg.iloc[i]['GENES']],
                    '# SNPS': [snp_list_ch2.iloc[j]['POS']]
                })
                disrupted_promoter_df = pd.concat([disrupted_promoter_df, new_entry], ignore_index=True)
    disrupted_promoter_df['# SNPS'] = disrupted_promoter_df.groupby(['OPERON_ID'])['# SNPS'].transform('count')
    disrupted_promoter_df = disrupted_promoter_df.drop_duplicates().reset_index(drop = True)
    
    disrupted_promoter_snp_df_list.append(disrupted_promoter_df)
    disrupted_promoter_snp_df_all = pd.concat([disrupted_promoter_snp_df_all, disrupted_promoter_df], ignore_index=True)

    n_operon = len(df_ch1) + len(df_ch2)
    n_snp = len(snp_list)
    n_broken_promoter = len(disrupted_promoter_df)
    p_broken_promoter = round(n_broken_promoter * 100 / n_operon, 1)
    broken_promoter = str(n_broken_promoter) + " ("+ str(p_broken_promoter) + "%)"
    
    disrupted_promoter_df = disrupted_promoter_df.assign(GENE=disrupted_promoter_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
    n_genes = len(disrupted_promoter_df)
    p_genes = round(n_genes * 100 / total_genes, 1)
    broken_genes = str(n_genes) + " ("+ str(p_genes) + "%)"    
    
    new_summary = pd.DataFrame.from_dict({
                        'SAMPLE': [sample],
                        '# SNPS': [n_snp],
                        'DISRUPTED PROMOTERS': [broken_promoter],
                        'DISRUPTED GENES': [broken_genes]
    })
    summary_promoter_snp = pd.concat([summary_promoter_snp, new_summary], ignore_index=True)


# Filter promoters to get Zagreb specific results
zagreb_promoter_snp_df = disrupted_promoter_snp_df_all[disrupted_promoter_snp_df_all['SAMPLE'].str.contains('Zagreb', na=False)].reset_index(drop = True)
zagreb_promoter_snp_df = zagreb_promoter_snp_df.assign(GENE=zagreb_promoter_snp_df.GENES.str.split(";")).explode('GENE').reset_index(drop = True)
zagreb_promoter_snp_df = pd.merge(zagreb_promoter_snp_df, df_rnaseq_final[['Locus_tag.gene', 'log2FC', 'qval']], left_on='GENE', right_on='Locus_tag.gene', how='left').reset_index(drop = True)
zagreb_promoter_snp_df = zagreb_promoter_snp_df.drop(['GENE', 'Locus_tag.gene'], axis=1)
zagreb_promoter_snp_df['log2FC_str'] = round(zagreb_promoter_snp_df['log2FC'],3).astype(str)
zagreb_promoter_snp_df['LOG2FC'] = zagreb_promoter_snp_df.groupby(['OPERON_ID'])['log2FC_str'].transform(lambda x: ';'.join(x))
zagreb_promoter_snp_df = zagreb_promoter_snp_df[zagreb_promoter_snp_df['qval'] < 1E-3]
zagreb_promoter_snp_df = zagreb_promoter_snp_df[abs(zagreb_promoter_snp_df['log2FC']) > 1]
zagreb_promoter_snp_df = zagreb_promoter_snp_df.drop(['log2FC', 'log2FC_str', 'qval'], axis=1).drop_duplicates().reset_index(drop = True)

# Write results to as xlsx sheets and csv files
writer = pd.ExcelWriter('promoter_snp_sheets_V3.xlsx', engine='xlsxwriter')
for i in range(len(disrupted_promoter_df_list)):
    disrupted_promoter_snp_df_list[i].to_excel(writer, sheet_name=sample_list[i], index=False)
writer.close()

summary_promoter_snp.to_csv('promoter_snp_summary_V3.csv', index=False)
disrupted_promoter_snp_df_all.to_csv('promoter_snp_all_V3.csv', index=False)
zagreb_promoter_snp_df.to_csv('zagreb_promoter_snp_rna_V3.csv', index=False)