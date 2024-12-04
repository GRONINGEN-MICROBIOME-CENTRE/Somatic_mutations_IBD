# load libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
from matplotlib.pyplot import figure
import itertools

# Read in vcf file output of mutect2 in a dataframe containing genes as rownames, biopsy samplename as columns
# and the number of somatic mutations as values
# this functionality is not very fast
#covered_genenames_wes.txt

# define all functions used for reading in files
def get_gene(values):
    # this function takes the info column and pulls out the gene name from the funcotator annotation
    return values.split('[')[1].split(']')[0].split('|')[0]


def pull_mutations(mutation_folder, phenotype_file, skipnrows=0):
    result_input = subprocess.run(['ls', mutation_folder], stdout=subprocess.PIPE)
    ls_input = result_input.stdout.decode("utf-8").split('\n')
    result_input = subprocess.run(['ls', "/Users/iwan/Research/Somatic_mutations/bedgraph_local"], stdout=subprocess.PIPE)
    read_depth_ls_input = result_input.stdout.decode("utf-8").split('\n')
    phenotype_df = pd.read_excel(new_phenotype_file, sheet_name=0).drop_duplicates(subset=['biopsy_number']).drop(["Unnamed: 42"], axis=1)
    #print(phenotype_df)
    biopsy_list_01 = phenotype_df['biopsy_number']
    year = phenotype_df['Sequencing_batch']
    # set index based on all biopsies
    index = biopsy_list_01
    #create one column with a gene name for creating an empty dataframe
    gene_name_list = pd.read_csv("/Users/iwan/Research/Somatic_mutations/25_percent_pipeline_output/gene_list_uniq_names.txt", header=None)[0].tolist()
    columns = gene_name_list
    # create and fill dataframe
    df_ = pd.DataFrame(index=index, columns=columns)
    #df_ = df_.fillna(0)
    for biopsy in biopsy_list_01:
    # find file matching biopsy from list of files
    # adapt this to also pull the coverage file with gene names of passes, then turn those values to zero.
        file = [file for file in ls_input if biopsy + '_' in file]
        file = [x for x in file if x.endswith('.vcf')]
        if file:
            #print(file)
            if file[0].startswith('control'):
                new_file = file
                #[new_file for new_file in file if biopsy in new_file.split('_')[1]]
            else:
                new_file = [new_file for new_file in file if biopsy in new_file.split('-')[1]]
            #print(new_file)
        # load file into pandas
            if new_file:
                sample_df = pd.read_csv('{}{}'.format(mutation_folder, new_file[0]), delimiter='\t', skiprows=skipnrows)
                sample_df['Genes'] = sample_df.iloc[:, 7].apply(get_gene)
        # check if all genes are already listed
                new_cols = list(set(sample_df['Genes'].tolist()) - set(df_.columns.tolist()))
                if new_cols:
                #this line makes pandas very unhappy, not sure how to fix
                    #df_[new_cols] = 0
                    df_[new_cols] = np.nan
                gene_numbers = sample_df['Genes'].value_counts()
                file = [file for file in read_depth_ls_input if biopsy + '_' in file]
                file = [x for x in file if x.endswith('_passed.txt')]
                # here we check if we have coverage over the gene, if yes we know there is no mutation
                if file:
                    gene_file = pd.read_csv("/Users/iwan/Research/Somatic_mutations/bedgraph_local/{}".format(file[0]), sep='\t', header=None)[2]
                    for gene in gene_file:
                        df_.at[biopsy, gene] = 0
                for gene in set(sample_df['Genes'].tolist()):
                    df_.at[biopsy, gene] = gene_numbers[gene]

    return df_

# weird inconsistency in the phenotype data

def BMI_fixer(x):
    new_bmi = [str(i).replace('.', '') for i in [x]]
    #print(new_bmi)
    new_bmi_2 = [str(i).replace(',', '') for i in new_bmi]
    newest_bmi = []
    #print(new_bmi_2)
    for value in new_bmi_2:
        if value != 'NA':
            if value != 'nan':
                if value != '':
                    newest_bmi.append(float('{}.{}'.format(value[:2], value[2:])))
                else:
                    newest_bmi.append('NA')
            else:
                    newest_bmi.append('NA')
        else:
                    newest_bmi.append('NA')
    return newest_bmi[0]

# homogenize phenotypic information
def sex_fixer(x):
    if x == "Female" or x == "female" or x == "Female ":
        return 2
    elif x == "Male" or x == "male" or x == "Male ":
        return 1
    else:
        return x



def Diagnosis_fixer(x):
    if x == "Control" or x == "control":
        return "Control"
    elif x == "CD" or x == "cd":
        return "CD"
    elif x == "UC" or x == "uc":
        return "UC"
    else:
        return x



new_phenotype_file = "/Users/iwan/Research/Somatic_mutations/metadata/Somatic_mutations_metadata.xlsx"

fully_filtered_loc = "/Users/iwan/Research/Somatic_mutations/all_vcf_output/raw_vcfs/"
df_somatic_mutations = pull_mutations(fully_filtered_loc, new_phenotype_file)
df_somatic_mutations_columns = df_somatic_mutations.columns

updated_somatic_columns = df_somatic_mutations.columns

phenotype_df = pd.read_excel(new_phenotype_file, sheet_name=0).drop_duplicates(subset=['biopsy_number']).drop(
    ["Unnamed: 42"], axis=1)
phenotype_df.index = phenotype_df['biopsy_number']


# here we fix all the phenotype data and add it to the dataframe

def montreal_L_fix(value):
    if type(value) == type(1):
        if int(value) == 0:
            return 1
        elif int(value) == 1:
            return 2
        elif int(value) == 2:
            return 3
        elif int(value) == 4:
            return 4
        elif int(value) == 5:
            return 5
        elif int(value) == 6:
            return 6
        else:
            return np.nan
    else:
        return np.nan


def montreal_B_fix(value):
    if type(value) == type(1):
        if int(value) == 0:
            return 1
        elif int(value) == 1:
            return 2
        elif int(value) == 2:
            return 3
        elif int(value) == 3:
            return 4
        elif int(value) == 4:
            return 5
        elif int(value) == 5:
            return 6
        else:
            return np.nan
    else:
        return np.nan


def montreal_A_fix(value):
    if value == 3:
        return 2
    elif value == 'A3':
        return 2
    elif value == 'A2':
        return 1
    elif value == 'A1':
        return 0
    else:
        return value


def montreal_E_fix(value):
    if value == 'E3':
        return 2
    elif value == 'E2':
        return 1
    elif value == 'E1':
        return 0
    else:
        return value


def montreal_S_fix(value):
    if value == 'S3':
        return 3
    elif value == 'S2':
        return 2
    elif value == 'S1':
        return 1
    elif value == 'S0':
        return 0
    else:
        return value

tnf_pre_list = []
tnf_during_list = []
tnf_combined = []
for med1, med2, med3 in zip(phenotype_df['IFX_use_all'], phenotype_df['ADA_use_all'], phenotype_df['TNF_other_all']):
    if (med1 == 1) or (med2 == 1) or (med3 == 1):
        tnf_pre_list.append(1)
        tnf_during_list.append(0)
        tnf_combined.append(1)
    elif (med1 == 2) or (med2 == 2) or (med3 == 2):
        tnf_during_list.append(1)
        tnf_pre_list.append(0)
        tnf_combined.append(1)
    elif (med1 == 3) or (med2 == 3) or (med3 == 3):
        tnf_pre_list.append(0)
        tnf_during_list.append(0)
        tnf_combined.append(0)
    elif (med1 == 4) or (med2 == 4) or (med3 == 4):
        tnf_pre_list.append(1)
        tnf_during_list.append(1)
        tnf_combined.append(1)
    elif (med1 == 5) or (med2 == 5) or (med3 == 5):
        tnf_during_list.append(1)
        tnf_pre_list.append(0)
        tnf_combined.append(1)
    elif (med1 == 6) or (med2 == 6) or (med3 == 6):
        tnf_pre_list.append(1)
        tnf_during_list.append(0)
        tnf_combined.append(1)
    elif (med1 == 7) or (med2 == 7) or (med3 == 7):
        tnf_pre_list.append(1)
        tnf_during_list.append(1)
        tnf_combined.append(1)
    else:
        tnf_pre_list.append(0)
        tnf_during_list.append(0)
        tnf_combined.append(0)

df_somatic_mutations['Montreal_L'] = phenotype_df['MontrealL']  # .apply(montreal_L_fix)
df_somatic_mutations['Montreal_B'] = phenotype_df['MontrealB']  # .apply(montreal_B_fix)
df_somatic_mutations['Montreal_E'] = phenotype_df['MontrealE']  # .apply(montreal_E_fix)
df_somatic_mutations['Montreal_S'] = phenotype_df['MontrealS']  # .apply(montreal_S_fix)
df_somatic_mutations['Montreal_A'] = phenotype_df['MontrealA']  # .apply(montreal_A_fix)
df_somatic_mutations['Diagnosis'] = phenotype_df['Diagnosis']
df_somatic_mutations['Sequencing_Year'] = phenotype_df['Sequencing_batch']
df_somatic_mutations['Inflammation'] = phenotype_df['Inflammation']
df_somatic_mutations['Combined'] = phenotype_df['Diagnosis'] + phenotype_df['Inflammation']
df_somatic_mutations['Age_at_biopsy'] = phenotype_df['age_at_biopsy']
df_somatic_mutations['Location_rough'] = phenotype_df['Location_rough']
df_somatic_mutations['BMI'] = phenotype_df['BMI']
df_somatic_mutations['Sex'] = phenotype_df['sex']
df_somatic_mutations['umcg_id'] = phenotype_df['Research ID']
df_somatic_mutations['tnf_during'] = tnf_during_list
df_somatic_mutations['smoking'] = phenotype_df['smoking_DB'].replace('Yes', 1).replace('No', 0).replace('yes',
                                                                                                        1).replace('no',
                                                                                                                   0)
df_somatic_mutations['Steroids'] = phenotype_df['Steroids'].replace('Yes', 1).replace('No', 0).replace('yes',
                                                                                                       1).replace('no',
                                                                                                                  0)
df_somatic_mutations['Thiopurines'] = phenotype_df['Thiopurines'].replace('Yes', 1).replace('No', 0).replace('yes',
                                                                                                             1).replace(
    'no', 0)
df_somatic_mutations['Methotrexaat'] = phenotype_df['Methotrexaat'].replace('Yes', 1).replace('No', 0).replace('yes',
                                                                                                               1).replace(
    'no', 0)
df_somatic_mutations['Aminosalicylates'] = phenotype_df['Aminosalicylates'].replace('Yes', 1).replace('No', 0).replace(
    'yes', 1).replace('no', 0)
# remove any rows with no mutations detected at all.


# add montreal A E S

df_somatic_mutations['BMI'] = df_somatic_mutations['BMI'].apply(BMI_fixer)
df_somatic_mutations['Diagnosis'] = df_somatic_mutations['Diagnosis'].apply(Diagnosis_fixer)
df_somatic_mutations['Sex'] = df_somatic_mutations['Sex'].apply(sex_fixer)

# df_rna_edit_filter = df_rna_edit_filter[(df_rna_edit_filter[df_rna_edit_filter_cols].T >= 0).any()]
df_somatic_mutations = df_somatic_mutations[(df_somatic_mutations[updated_somatic_columns].fillna(0).T > 0).any()]

# store output
df_somatic_mutations = df_somatic_mutations.replace('Light', 'Yes')
df_somatic_mutations.to_csv('/Users/iwan/Research/Somatic_mutations/output_data/Somatic_mutations_with_controls.csv')

# this is repeated for the other input datasets, e.g. containing only pathogenic mutations




# we correct for gene length in one of the analysis:
# this code is used for this correction
def genelength_adder(gene):
    result_input = subprocess.run(["grep", "{}".format(gene), "/Users/iwan/Research/Somatic_mutations/metadata/Homo_sapiens.GRCh38.104.chr.gff3"], stdout=subprocess.PIPE).stdout.decode("utf-8").split('\n')
    result_input = [x for x in result_input if 'ensembl_havana\tgene' in x]
    if result_input == []:
        return np.nan
    else:
        genelength = int(result_input[0].split('\t')[4]) - int(result_input[0].split('\t')[3])
        return genelength

gene_length_dict = {}


for genename in df_somatic_mutations[updated_somatic_columns].columns:
    gene_length_dict[genename] = genelength_adder(genename)

df_somatic_corrected = df_somatic_mutations.copy()
for gene in updated_somatic_columns:
    df_somatic_corrected[gene] = df_somatic_corrected[gene] / gene_length_dict[gene]


