import sys
import pandas as pd
import numpy as np
import pysam
import gzip
from itertools import groupby

argv = sys.argv[1:]
path_bam = argv[0]  # input bam
path_vcf= argv[1]
basename = argv[2]      # output path
chromosome = argv[3]
gene_start = int(argv[4])
gene_stop = int(argv[5])
min_variants_phase_block = int(argv[6])-1
viewpointpos = int(argv[7])
fill_undifinded_varaints = bool(argv[8]=='True')

# test

#path_bam =  '../phasing/nanopore/CYP2C19_chr10_filter_100_phased_2.bam'
#path_vcf = '../phasing/nanopore/CYP2C19_chr10_phased_haplo.vcf.gz'
#basename =  '../phasing/nanopore/CYP2C19_2_'
#chromosome = 'chr10'
#gene_start = 94761287
#gene_stop = 94853205
#min_variants_phase_block = 1
#viewpointpos = 0
#fill_undifinded_varaints = True

print('Start haplotyping')
print('Input:')
print(f"  Input bam:       {path_bam}")
print(f"  Output path:     {basename}")
print(f"  Chromosome:      {chromosome}")
print(f"  Gene start:      {gene_start}")
print(f"  Gene stop:       {gene_stop}")
print(f"  Min variants:    {min_variants_phase_block}")
print(f"  Adjust variants: {fill_undifinded_varaints}")

def recalculate_all(hap1, hap2):
    allpositions = pd.concat([hap1['pos'], hap2['pos']], axis=0).drop_duplicates()
    for position in allpositions:
        searchID1 = hap1[hap1['pos']==position]['ID']
        searchID2 = hap2[hap2['pos']==position]['ID']
        
        #check for the position with the searchID how strong it's linked to the current haplotype 1
        if len(hap1[hap1['pos']==position])==1:
            #extract all lines with where ID1 matches the searchID
            links = hapfile[hapfile['ID1'].isin(searchID1)]
            #extract all of these where the other ID is in either haplotype (split)
            linkstohap1 = links[links['ID2'].isin(hap1['ID'])]
            linkstohap2 = links[links['ID2'].isin(hap2['ID'])]

            #get according weight/plexity values
            weight_1 = linkstohap1['weight'].sum()
            ambig_1 = linkstohap2['weight'].sum()
            plexity_1 = len(linkstohap1)
            amplexity_1 = len(linkstohap2)
        else:
            weight_1 = 0
            ambig_1 = 0
            plexity_1 = 0
            amplexity_1 = 0
        
        #check for the position with the searchID how strong it's linked to the current haplotype 2
        if len(hap2[hap2['pos']==position]==1):
            #extract all lines with where ID1 matches the searchID
            links = hapfile[hapfile['ID1'].isin(searchID2)]
            #extract all of these where the other ID is in either haplotype (split)
            linkstohap2 = links[links['ID2'].isin(hap2['ID'])]
            linkstohap1 = links[links['ID2'].isin(hap1['ID'])]
            
            #get according weight/plexity values
            weight_2 = linkstohap2['weight'].sum()
            ambig_2 = linkstohap1['weight'].sum()
            plexity_2 = len(linkstohap2)
            amplexity_2 = len(linkstohap1)
        else:
            weight_2 = 0
            ambig_2 = 0
            plexity_2 = 0
            amplexity_2 = 0
        
        #get final recalculation values for the current snp
        weight = weight_1 + weight_2
        ambig = ambig_1 + ambig_2
        plexity = plexity_1 + plexity_2
        amplexity = amplexity_1 + amplexity_2
        ratio = round(ambig/weight, 3)

        #assign recalculated values depending on hap1 or hap2
        hap1.loc[hap1.pos == position, 'weight'] = weight
        hap1.loc[hap1.pos == position, 'ambig'] = ambig
        hap1.loc[hap1.pos == position, 'plexity'] = plexity
        hap1.loc[hap1.pos == position, 'amplexity'] = amplexity
        hap1.loc[hap1.pos == position, 'ratio'] = ratio
        if (weight_1 > 0 and weight_2 > 0 and(hap1.loc[hap1.pos == position, 'class'].values[0] != "VP")):
            hap1.loc[hap1.pos == position, 'class'] = "Double_link"

        hap2.loc[hap2.pos == position, 'weight'] = weight
        hap2.loc[hap2.pos == position, 'ambig'] = ambig
        hap2.loc[hap2.pos == position, 'plexity'] = plexity
        hap2.loc[hap2.pos == position, 'amplexity'] = amplexity
        hap2.loc[hap2.pos == position, 'ratio'] = ratio
        if (weight_1 > 0 and weight_2 > 0 and(hap2.loc[hap2.pos == position, 'class'].values[0] != "VP")):
            hap2.loc[hap2.pos == position, 'class'] = "Double_link"

    return(hap1, hap2)
            
def target_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    read_consuming_ops = ("M", "D", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
    return result


#Create names
hap1_out_file = f"{basename}{chromosome}_hap1.txt"
hap2_out_file = f"{basename}{chromosome}_hap2.txt"
haps_out_file = f"{basename}{chromosome}_haps.txt"
netwtable_file = f"{basename}{chromosome}_network.txt"
vcf_out = f"{basename}{chromosome}_phased_haplo.vcf"

print("Make link table")


file = gzip.open(path_vcf,'rt')
header_count = 0


with open(vcf_out,'w') as fh:
    for line in file.readlines():
        if line.startswith("##"):
            header_count += 1
            fh.write(line)
        else:
            break

file.close()
variants = pd.read_csv(path_vcf,sep='\t',header=header_count)
variants.index = variants['POS']
selected_variants = variants.loc[gene_start:gene_stop]

haplo_part = 0
haplo_dict = dict()
haplo_count = 0
haplo_pos = []
for row in variants.index:
    split_sample = variants.loc[row,'SAMPLE'].split(':')
    if len(split_sample) == 5 and ('|' in split_sample[0]):
        if haplo_part == split_sample[4]:
            haplo_count += 1
            haplo_pos.append(int(row))
        else:
            if haplo_count > min_variants_phase_block:
                haplo_dict[int(haplo_part)]={'pos':haplo_pos}
            haplo_part = split_sample[4]
            haplo_count = 1
            haplo_pos = [int(row)]

for variant in haplo_dict:
    haplo_dict[variant]['filter'] = False
    if len(haplo_dict[variant]['pos']) > 10:
        for i in range(len(haplo_dict[variant]['pos']) - 10):
            if haplo_dict[variant]['pos'][i] - haplo_dict[variant]['pos'][i+10] < 100:
                haplo_dict[variant]['filter'] = True

bam_file = pysam.AlignmentFile(path_bam)
read_dict = dict()

for read in bam_file.fetch():
    if read.reference_name == chromosome:
        if (gene_start < (read.reference_start + target_len(read.cigarstring))) and (read.reference_start < gene_stop):
            tags_dict = dict(read.tags)
            if 'HP' in tags_dict:
                if tags_dict['PS'] in haplo_dict:
                    if not haplo_dict[tags_dict['PS']]['filter']:
                        add = False
                        for i in range(len(haplo_dict[tags_dict['PS']]['pos'])-min_variants_phase_block):
                            start = haplo_dict[tags_dict['PS']]['pos'][i]
                            stop = haplo_dict[tags_dict['PS']]['pos'][i+min_variants_phase_block]
                            if (stop < (read.reference_start + target_len(read.cigarstring))) and (read.reference_start < start):
                                add = True
                        if add:
                            name = '_'.join(read.qname.split('_')[0:2])
                            if name in read_dict:
                                if tags_dict['PS'] in read_dict[name]:
                                    if  read_dict[name][tags_dict['PS']] != tags_dict['HP']:
                                        read_dict[name]['remove']=True
                                else:
                                    read_dict[name][tags_dict['PS']] = tags_dict['HP']
                            else:
                                read_dict[name] = {tags_dict['PS']:tags_dict['HP'],
                                                   'remove':False}

read_dict_clean = dict()
for read in read_dict:
    if not read_dict[read]['remove'] and (len(read_dict[read])>2):
        read_dict_clean[read] = read_dict[read]

if read_dict_clean != dict():
    
    linked_table = pd.DataFrame(columns = ['hap1','pos1','hap2','pos2','weigth'])
    
    for read in read_dict_clean:
        for pos1 in read_dict_clean[read]:
            for pos2 in read_dict_clean[read]:
                if (pos1 != pos2) and (pos1 != 'remove') and (pos2 != "remove"):
                    sample = '_'.join([str(read_dict_clean[read][pos1]),
                                           str(pos1),
                                           str(read_dict_clean[read][pos2]),
                                           str(pos2)])
                    if sample in linked_table.index:
                        linked_table.loc[sample,'weigth'] += 1
                    else:
                        linked_table.loc[sample] = [read_dict_clean[read][pos1],
                                                    pos1,
                                                    read_dict_clean[read][pos2],
                                                    pos2,1]
    
    bam_file.close()
    
    hapfile = linked_table.astype({'hap1':str,'hap2':str})
    hapfile.columns = ['var1', 'pos1', 'var2', 'pos2', 'weight']
    
    print("Start phasing")
    
    #If there's no viewpoint override, autodefine best VP based on how many reads with links were found for each position
    weight_frame = hapfile[['pos1','weight']]
    weight_frame['total']=weight_frame.groupby(['pos1'])['weight'].transform('sum')
    weight_frame = weight_frame[['pos1','total']].drop_duplicates()
    max_weight = weight_frame[weight_frame['total'] == weight_frame['total'].max()]
    
    if viewpointpos == 0:
        viewpointpos = int(max_weight['pos1'].head(1))
    
    viewpointvars='1','2'
    
    #create hap1 and assign the VP values
    hap1 = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
    
    hap1 = pd.concat([hap1,pd.DataFrame({
        'pos' : viewpointpos,
        'var' : viewpointvars[0],
        'ID' : (str(viewpointpos) + ":" + viewpointvars[0]),
        'weight' : 10000,
        'linkedweight' : 10000,
        'ambig' : 0,
        'class' : "VP",
        'ratio' : 0,
        'round' : 0,
        'plexity' : 0,
        'amplexity' : 0
        }, index=[0])])
    
    #create hap2 and assign the VP values
    hap2 = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
    hap2 = pd.concat([hap2,pd.DataFrame({
        'pos' : viewpointpos,
        'var' : viewpointvars[1],
        'ID' : (str(viewpointpos) + ":" + viewpointvars[1]),
        'weight' : 10000,
        'linkedweight' : 10000,
        'ambig' : 0,
        'class' : "VP",
        'ratio' : 0,
        'round' : 0,
        'plexity' : 0,
        'amplexity' : 0
        }, index=[0])])
    
    #initiate seedhap for first round
    seedhap1 = hap1
    seedhap2 = hap2
    
    #create IDlist with all positional values of the SNVs
    IDlist = hapfile['pos1'].drop_duplicates()
    IDlist = IDlist[IDlist != viewpointpos]
    IDlistnew = IDlist
    
    hapfile['ID1']=hapfile["pos1"].astype(str)+":"+hapfile['var1']
    hapfile['ID2']=hapfile["pos2"].astype(str)+":"+hapfile['var2']
    
    iterations = 25
    iteration = iter(range(iterations))
    average_ratio_1 = 0
    
    #Loop to build on haplotype, each time using the haplotypes as seeds
    for i in iteration:
        loopnr = i+1
    
    
        if IDlistnew.size > 0:
            badsnps = {"isolated": [], "single_allele": [], "ambigious": []}
            for snp in IDlistnew:
                #all lines from hapfile which contain the current SNP position were itterating over
                links = hapfile[hapfile['pos1']==snp]
                #all lines from that which have the other position in the current haplotype 1
                hap_links = links[links['pos2'].isin(seedhap1['pos'])]
                
                if (len(links['ID1'].unique())==2):
                    #get all links that belong to one or the other ID for the linked variant
                    var1 = links['ID1'].drop_duplicates().reset_index(drop=True)[0]
                    var2 = links['ID1'].drop_duplicates().reset_index(drop=True)[1]
                    var1_links = hap_links[hap_links['ID1'] == var1]
                    var2_links = hap_links[hap_links['ID1'] == var2]
    
                    #determine weight for possibilities of how vars can be linked to haplotypes
                    #posibility 1: allele/hap 1/1 - 2/2
                    var1_hap1 = var1_links[var1_links["ID2"].isin(seedhap1['ID'])]
                    var2_hap2 = var2_links[var2_links["ID2"].isin(seedhap2['ID'])]
                    p1_weight = var1_hap1['weight'].sum() + var2_hap2['weight'].sum()
                    p1_plexity = len(var1_hap1) + len(var2_hap2)
    
                    #posibility 1: allele/hap 1/2 - 2/1
                    var1_hap2 = var1_links[var1_links["ID2"].isin(seedhap2['ID'])]
                    var2_hap1 = var2_links[var2_links["ID2"].isin(seedhap1['ID'])]
                    p2_weight = var1_hap2['weight'].sum() + var2_hap1['weight'].sum()
                    p2_plexity = len(var1_hap2) + len(var2_hap1)
    
                    if(p1_weight == 0 and p2_weight == 0):
                        badsnps["isolated"].append(snp)
    
                    if(p1_weight == p2_weight):
                        badsnps["ambigious"].append(snp)
                    
                    #if var1 is in hap1 en vice versa
                    if(p1_weight>p2_weight):
                        weight = p1_weight
                        ambig = p2_weight
                        plexity = p1_plexity
                        amplexity = p2_plexity
                    
                        #check if it's a double of single link
                        if (var1_hap1['weight'].sum() > 0 and var2_hap2['weight'].sum() > 0):
                            snv_class_1 = "Double_link"
                            snv_class_2 = "Double_link"
                        elif (var1_hap1['weight'].sum() > 0 and var2_hap2['weight'].sum() == 0): 
                            snv_class_1 = "Single_link"
                            snv_class_2 = "Indirect"
                        elif (var2_hap2['weight'].sum() > 0 and var1_hap1['weight'].sum() == 0): 
                            snv_class_2 = "Single_link"
                            snv_class_1 = "Indirect"
    
                        #write line for haplotype 1
                        hap1 = pd.concat([hap1,pd.DataFrame({
                                    'pos' : snp,
                                    'var' : var1[-1:],
                                    'ID' : var1,
                                    'weight' : weight,
                                    'linkedweight' : weight,
                                    'ambig' : ambig,
                                    'class' : snv_class_1,
                                    'ratio' : round(ambig/weight, 3),
                                    'round' : i+1,
                                    'plexity' : plexity,
                                    'amplexity' : amplexity
                                    }, index=[0])])
    
                        #write line for haplotype 2
                        hap2 = pd.concat([hap2,pd.DataFrame({
                                    'pos' : snp,
                                    'var' : var2[-1:],
                                    'ID' : var2,
                                    'weight' : weight,
                                    'linkedweight' : weight,
                                    'ambig' : ambig,
                                    'class' : snv_class_2,
                                    'ratio' : round(ambig/weight, 3),
                                    'round' : i+1,
                                    'plexity' : plexity,
                                    'amplexity' : amplexity
                                    }, index=[0])])
    
                    #if var1 is in hap2 en vice versa
                    if(p2_weight>p1_weight):
                        weight = p2_weight
                        ambig = p1_weight
                        plexity = p2_plexity
                        amplexity = p1_plexity
                    
                        #check if it's a double of single link
                        if (var1_hap2['weight'].sum() > 0 and var2_hap1['weight'].sum() > 0):
                            snv_class_1 = "Double_link"
                            snv_class_2 = "Double_link"
                        elif (var2_hap1['weight'].sum() > 0 and var1_hap2['weight'].sum() == 0): 
                            snv_class_1 = "Single_link"
                            snv_class_2 = "Indirect"
                        elif (var1_hap2['weight'].sum() > 0 and var2_hap1['weight'].sum() == 0): 
                            snv_class_2 = "Single_link"
                            snv_class_1 = "Indirect"
    
                        #write line for haplotype 1
                        hap1 = pd.concat([hap1,pd.DataFrame({
                                    'pos' : snp,
                                    'var' : var2[-1:],
                                    'ID' : var2,
                                    'weight' : weight,
                                    'linkedweight' : weight,
                                    'ambig' : ambig,
                                    'class' : snv_class_1,
                                    'ratio' : round(ambig/weight,3),
                                    'round' : loopnr,
                                    'plexity' : plexity,
                                    'amplexity' : amplexity
                                    }, index=[0])])
    
                        #write line for haplotype 2
                        hap2 = pd.concat([hap2,pd.DataFrame({
                                    'pos' : snp,
                                    'var' : var1[-1:],
                                    'ID' : var1,
                                    'weight' : weight,
                                    'linkedweight' : weight,
                                    'ambig' : ambig,
                                    'class' : snv_class_2,
                                    'ratio' : round(ambig/weight,3),
                                    'round' : loopnr,
                                    'plexity' : plexity,
                                    'amplexity' : amplexity
                                    }, index=[0])])
    
                else:
                    badsnps["single_allele"].append(snp)
        
        #Subselect and recalculate with decreasing stringency per iteration
        if loopnr < 5:
            quantileweight = 30-(5*loopnr)
            hap1 = hap1[(hap1['ratio']<0.1)|(hap1['class']=="VP")]
            hap2 = hap2[(hap2['ratio']<0.1)|(hap2['class']=="VP")]
    
            hap1, hap2 = recalculate_all(hap1, hap2)
    
            hap1 = hap1[(hap1['linkedweight']>quantileweight)&(hap1['class']!="Single_link")&(hap1['class']!="Indirect")|(hap1['class']=="VP")]
            hap2 = hap2[(hap2['linkedweight']>quantileweight)&(hap2['class']!="Single_link")&(hap2['class']!="Indirect")|(hap2['class']=="VP")]
    
        elif loopnr < (iterations-4): 
            hap1 = hap1[hap1['ratio']<0.3]
            hap2 = hap2[hap2['ratio']<0.3]
            
            hap1, hap2 = recalculate_all(hap1, hap2)
    
            hap1 = hap1[(hap1['class']!="Single_link")&(hap1['class']!="Indirect")]
            hap2 = hap2[(hap2['class']!="Single_link")&(hap2['class']!="Indirect")]
            
        else:
            hap1, hap2 = recalculate_all(hap1, hap2)
    
        #Update the seedhap and IDlist for the next round
        seedhap1 = hap1
        seedhap2 = hap2
    
        IDlistnew = IDlist[~IDlist.isin(seedhap1['pos'])]
        IDlistnew = IDlistnew[~IDlistnew.isin(seedhap2['pos'])]
    
    #Re-evaluate bad snps
    final_badsnps = {"Single_linked": [], "Single_ambiguous": [], "Single_isolated": [], "Isolated": []}
    
    hap1_b = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
    hap2_b = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
    
    for snp in badsnps["single_allele"]:
        #get all links that with the snp that link to hap1
        links = hapfile[hapfile['pos1'] == snp]
        hap_links = links[links['pos2'].isin(hap1['pos'])]
        
        #Split all these links into whether they would fit in hap1 or hap2
        hap1_links = hap_links[hap_links['ID2'].isin(hap1['ID'])]
        hap2_links = hap_links[hap_links['ID2'].isin(hap2['ID'])]
    
        #get the NT for each to use when adding to the dataframe
        hap1_var = hap1_links['var1'].drop_duplicates().reset_index(drop=True)
        hap2_var = hap2_links['var1'].drop_duplicates().reset_index(drop=True)
    
        #get weight and plexity for each    
        weight_1 = hap1_links['weight'].sum()
        weight_2 = hap2_links['weight'].sum()        
        plexity_1 = len(hap1_links)
        plexity_2 = len(hap2_links)
    
        #determine whether the badsnp fits into one, both or none of the haplotypes and process accordingly
        if weight_1 > 0 and weight_2 == 0:
            final_badsnps["Single_linked"].append(snp)
            hap1 = pd.concat([hap1,pd.DataFrame({
                        'pos' : snp,
                        'var' : hap1_var[0],
                        'ID' : str(snp) + ":" + hap1_var[0],
                        'weight' : weight_1,
                        'linkedweight' : weight_1,
                        'ambig' : weight_2,
                        'class' : "Single_link",
                        'ratio' : round(weight_2/weight_1,3),
                        'round' : 26,
                        'plexity' : plexity_1,
                        'amplexity' : plexity_2
                        }, index=[0])])
    
        if weight_2 > 0 and weight_1 == 0:
            final_badsnps["Single_linked"].append(snp)
            hap2 = pd.concat([hap2,pd.DataFrame({
                        'pos' : snp,
                        'var' : hap2_var[0],
                        'ID' : str(snp) + ":" + hap2_var[0],
                        'weight' : weight_2,
                        'linkedweight' : weight_2,
                        'ambig' : weight_1,
                        'class' : "Single_link",
                        'ratio' : round(weight_1/weight_2,3),
                        'round' : 26,
                        'plexity' : plexity_2,
                        'amplexity' : plexity_1
                        }, index=[0])])
                        
        if weight_1 > 0 and weight_2 > 0:
            final_badsnps["Single_ambiguous"].append(snp)
    
        if weight_1 == 0 and weight_2 == 0:
            final_badsnps["Single_isolated"].append(snp)
    
    #Filter bad snps for acceptable ones after recalculation
    hap1_c = pd.concat([hap1, hap1_b], axis=0)    
    hap2_c = pd.concat([hap2, hap2_b], axis=0)    
    
    hap1_c, hap2_c = recalculate_all(hap1=hap1_c, hap2=hap2_c)
    
    hap1_c = hap1_c[hap1_c['ratio']<0.5]
    hap2_c = hap2_c[hap2_c['ratio']<0.5]
    
    hap1_b = hap1_b[hap1_b['ID'].isin(hap1_c['ID'])]
    hap2_b = hap2_b[hap2_b['ID'].isin(hap2_c['ID'])]
    
    hap1 = pd.concat([hap1, hap1_b], axis=0)
    hap2 = pd.concat([hap2, hap2_b], axis=0)
    
    #Final recalculation and filtering for acceptable ratio
    hap1, hap2 = recalculate_all(hap1, hap2)
    
    #get all single_linked_bad positions
    hap1_sl_pos = hap1[(hap1['class']=="Single_link")&(~hap1['pos'].isin(hap2['pos']))]['pos']
    hap2_sl_pos = hap2[(hap2['class']=="Single_link")&(~hap2['pos'].isin(hap1['pos']))]['pos']
    
    for pos in hap1_sl_pos:
        #get both possible vars at position from pgsnp file
        var1 = hap1[hap1['pos']==pos]['var'].values[0]
        
        #find which is var for each
        if '1' == var1:
            var_hap2 = '2'
        if '2' == var1:
            var_hap2 = '1'
    
        #copy hap2 line and adjust to make hap1 indirect linked snp
        hap2_line = hap1[hap1['pos']==pos]
        hap2_line['var'] = var_hap2
        hap2_line['ID'] = str(pos) + ":" + var_hap2
        hap2_line['class'] = "Indirect"
    
        #add to hap1
        hap2 = pd.concat([hap2, hap2_line])
    
    for pos in hap2_sl_pos:
        #get both possible vars at position from pgsnp file
        var2 = hap2[hap2['pos']==pos]['var'].values[0]
        
        #find which is var for each
        if '1' == var2:
            var_hap1 = '2'
        if '2' == var2:
            var_hap1 = '2'
    
        #copy hap2 line and adjust to make hap1 indirect linked snp
        hap1_line = hap2[hap2['pos']==pos]
        hap1_line['var'] = var_hap1
        hap1_line['ID'] = str(pos) + ":" + var_hap1
        hap1_line['class'] = "Indirect"
    
        #add to hap1
        hap1 = pd.concat([hap1, hap1_line])
    
    hap1 = hap1[hap1['ratio']<0.5]
    hap2 = hap2[hap2['ratio']<0.5]
    
    hap1.drop_duplicates(subset=['pos'], inplace=True, keep='last')
    hap2.drop_duplicates(subset=['pos'], inplace=True, keep='last')
    
    hap1['hap'] = 1
    hap2['hap'] = 2
    
    #Write output files
    hap1_out = hap1[['pos', 'var']]
    hap2_out = hap2[['pos', 'var']]
    haps_out = pd.concat([hap1, hap2], axis=0)
    
    hap1_out.to_csv(hap1_out_file, sep=":", index=False, header=False)
    hap2_out.to_csv(hap2_out_file, sep=":", index=False, header=False)
    haps_out.to_csv(haps_out_file, sep="\t", index=False)
    
    #make network files
    
    hapfile['node1added'] = "X"
    hapfile['node2added'] = "X"
    
    hap1_links_nw = ""
    hap2_links_nw = ""
    
    for index, row in hap1.iterrows():
        connectables = hap1[hap1['round']<row['round']]['ID']
        links_nw = hapfile[hapfile['ID2']==row['ID']]
        links_nw = links_nw[links_nw['ID1'].isin(connectables)]
    
        if len(links_nw)>0:
            links_nw['used'] = row['round']
            links_nw['node2added'] = row['round']
            for index, row in links_nw.iterrows():
                
                linkID = links_nw.loc[[index],['ID1']].values[0][0]
                if len(hap1[hap1['ID']==linkID]['round']) > 0:
                    node1added = hap1[hap1['ID']==linkID]['round'].values[0]
                    links_nw.loc[[index],['node1added']] = node1added
    
            if len(hap1_links_nw) > 0:
                hap1_links_nw = pd.concat([hap1_links_nw, links_nw], axis=0)
            else:
                hap1_links_nw = links_nw  
        else:
            hap1.loc[[index],['class']] = "indirect"
    
    for index, row in hap2.iterrows():
        connectables = hap2[hap2['round']<row['round']]['ID']
        links_nw = hapfile[hapfile['ID1']==row['ID']]
        links_nw = links_nw[links_nw['ID2'].isin(connectables)]
    
        if len(links_nw)>0:
            links_nw['used'] = row['round']
            links_nw['node2added'] = row['round']
            for index, row in links_nw.iterrows():
                linkID = links_nw.loc[[index],['ID1']].values[0][0]
            
                if len(hap2[hap2['ID']==linkID]['round']) > 0:
                    node1added = hap2[hap2['ID']==linkID]['round'].values[0]
                    links_nw.loc[[index],['node1added']] = node1added
    
            if len(hap2_links_nw) > 0:
                hap2_links_nw = pd.concat([hap2_links_nw, links_nw], axis=0)
            else:
                hap2_links_nw = links_nw  
        else:
            hap2.loc[[index],['class']] = "indirect"
    
    hap1_links_nw['hap'] = 1
    hap1_links_nw['hap.1'] = 1
    hap1_links_nw['interaction'] = "pp"
    
    hap2_links_nw['hap'] = 2
    hap2_links_nw['hap.1'] = 2
    hap2_links_nw['interaction'] = "pp"
    
    all_links = pd.concat([hap1_links_nw, hap2_links_nw], axis = 0)
    
    both = pd.concat([hap1, hap2])
    
    
    netwtable = all_links.iloc[:, [5,6,12,9,7,8,10,11,1,3,4]]
    
    netwtable['weight1'] = "X"
    netwtable['weight2'] = "X"
    netwtable['ambig1'] = "X"
    netwtable['ambig2'] = "X"
    netwtable['plexity1'] = "X"
    netwtable['plexity2'] = "X"
    netwtable['amplexity1'] = "X"
    netwtable['amplexity2'] = "X"
    
    for index, row in netwtable.iterrows():
        hapdata1 = both[both['ID'] == row['ID1']].reset_index(drop=True)
        if len(hapdata1) == 1:
            netwtable.loc[[index],['weight1']] = hapdata1['weight'].values[0]
            netwtable.loc[[index],['ambig1']] = hapdata1['ambig'].values[0]
            netwtable.loc[[index],['plexity1']] = hapdata1['plexity'].values[0]
            netwtable.loc[[index],['amplexity1']] = hapdata1['amplexity'].values[0]
    
        
        hapdata2 = both[both['ID'] == row['ID2']].reset_index(drop=True)
        if len(hapdata2) == 1:
            netwtable.loc[[index],['weight2']] = hapdata2['weight'].values[0]
            netwtable.loc[[index],['ambig2']] = hapdata2['ambig'].values[0]
            netwtable.loc[[index],['plexity2']] = hapdata2['plexity'].values[0]
            netwtable.loc[[index],['amplexity2']] = hapdata2['amplexity'].values[0]
    
    
    netwtable['edgeweight'] = netwtable['weight']
    netwtable.drop(columns = ['weight'], inplace = True)
    
    netwtable.columns = ["ID1", "ID2", "interaction","interactionused", "added", "added", "hap", "hap",
                        "pos", "pos","weight", "weight", "ambig", "ambig", "plexity", "plexity", "amplexity", "amplexity", 
                        "interactionweight"]
    
    netwtable.to_csv(netwtable_file, sep="\t", index=False)
    
    big_haplo_group = list(hap1.loc[[x=='1' for x in hap1['var']]]['pos'])[0]
    
else:
    haplogroups = []
    for haplogroup in haplo_dict:
        if sum([1 for x in haplo_dict[haplogroup]['pos'] if gene_start < x < gene_stop]) > 0:
            haplogroups.append(haplogroup)
    if haplogroups == []:
        for row in  variants.loc[gene_start:gene_stop].index:
            if selected_variants.loc[row,'FORMAT'].split(':')[-1] == 'PS' and selected_variants.loc[row,'SAMPLE'].split(':')[-1] != '.':
                haplogroups.append(selected_variants.loc[row,'SAMPLE'].split(':')[-1])
        big_haplo_group = max(set(haplogroups),key=haplogroups.count)
    else:
        big_haplo_group = haplogroups[np.argmax([len(haplo_dict[haplogroup]['pos']) for haplogroup in haplogroups])]
    if haplogroups == []:
        big_haplo_group = 0
    hap1=pd.DataFrame({'pos':[big_haplo_group],'var':['1']})
    hap2=pd.DataFrame({'pos':[big_haplo_group],'var':['2']})
                                        
print("Split bam")

bam1_path = f"{basename}{chromosome}_hap1.bam"
bam2_path = f"{basename}{chromosome}_hap2.bam"

bam_file = pysam.AlignmentFile(path_bam)

hap1_reads=[]
hap2_reads=[]

for read in bam_file.fetch():
    if read.reference_name == chromosome:
        tags_dict = dict(read.tags)
        if 'HP' in tags_dict:
            if tags_dict['PS'] in list(hap1['pos']):
                if str(tags_dict['HP']) == list(hap1[hap1['pos'] == tags_dict['PS']]['var'])[0]:
                    hap1_reads.append('_'.join(read.qname.split('_')[0:2]))
            if tags_dict['PS'] in list(hap2['pos']):       
                if str(tags_dict['HP']) == list(hap2[hap2['pos'] == tags_dict['PS']]['var'])[0]:
                    hap2_reads.append('_'.join(read.qname.split('_')[0:2]))


hap1_set = set(hap1_reads) - set(hap2_reads)
hap2_set = set(hap2_reads) - set(hap1_reads)

bam_file.close()

bam_file = pysam.AlignmentFile(path_bam)

bam1_file = pysam.AlignmentFile(bam1_path, "wb", template=bam_file)
bam2_file = pysam.AlignmentFile(bam2_path, "wb", template=bam_file)

for read in bam_file.fetch():
    if '_'.join(read.qname.split('_')[0:2]) in hap1_set:
        bam1_file.write(read)
    if '_'.join(read.qname.split('_')[0:2]) in hap2_set:
        bam2_file.write(read)   
 
bam1_file.close()
bam2_file.close()
bam_file.close()

pysam.index(bam1_path)
pysam.index(bam2_path)


undefined_variants = []
for row in selected_variants.index:
    if selected_variants.loc[row,'FORMAT'].split(':')[-1] == 'PS':
        split_list = selected_variants.loc[row,'SAMPLE'].split(':')
        if split_list[-1] != '.':
            if int(split_list[-1]) in list(hap1['pos']):
                split_list[-1] = str(big_haplo_group)
                if list(hap1[hap1['pos'] == int(split_list[-1])]['var'])[0] != '1':
                    split_list[0] = '|'.join(split_list[0].split('|')[::-1])
                selected_variants.loc[row,'SAMPLE'] = ':'.join(split_list)
            else:
                undefined_variants.append(row)
        else:
            undefined_variants.append(row)
    else:
        undefined_variants.append(row)
            
if fill_undifinded_varaints:            
    
    if undefined_variants != []:
        bam1_file = pysam.AlignmentFile(bam1_path)
        
        phase_variant_dict = dict()
        for pileupcolumn in bam1_file.pileup():     
            if pileupcolumn.pos+1 in  undefined_variants:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:  
                        if pileupcolumn.pos+1 in phase_variant_dict:
                            if pileupread.alignment.query_sequence[pileupread.query_position] in phase_variant_dict[pileupcolumn.pos+1][1]:
                                phase_variant_dict[pileupcolumn.pos+1][1][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                            else:
                                phase_variant_dict[pileupcolumn.pos+1][1][pileupread.alignment.query_sequence[pileupread.query_position]] = 1
                        else:
                            phase_variant_dict[pileupcolumn.pos+1]={1:{pileupread.alignment.query_sequence[pileupread.query_position]:1}}
        
                
        bam2_file = pysam.AlignmentFile(bam2_path)
    
        for pileupcolumn in bam2_file.pileup():     
            if pileupcolumn.pos+1 in  undefined_variants:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:  
                        if pileupcolumn.pos+1 in phase_variant_dict:
                            if 2 in phase_variant_dict[pileupcolumn.pos+1]:
                                if pileupread.alignment.query_sequence[pileupread.query_position] in phase_variant_dict[pileupcolumn.pos+1][2]:
                                    phase_variant_dict[pileupcolumn.pos+1][2][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                                else:
                                    phase_variant_dict[pileupcolumn.pos+1][2][pileupread.alignment.query_sequence[pileupread.query_position]] = 1
                            else:
                                phase_variant_dict[pileupcolumn.pos+1][2] = {pileupread.alignment.query_sequence[pileupread.query_position]:1}
                        else:
                            phase_variant_dict[pileupcolumn.pos+1]={2:{pileupread.alignment.query_sequence[pileupread.query_position]:1}}      
                
        bam1_file.close()
        bam2_file.close()
        
    
        for pos in phase_variant_dict:
            split_format = selected_variants.loc[pos,'FORMAT'].split(':')
            split_list = selected_variants.loc[pos,'SAMPLE'].split(':')
            if split_format[-1] == 'PS':
                split_list[-1] = str(big_haplo_group)
            else:
                split_format.append('PS')
                split_list.append(str(big_haplo_group))
                selected_variants.loc[pos,'FORMAT'] = ':'.join(split_format)
            phase_haplo = []
            if 1 in phase_variant_dict[pos]:
                if selected_variants.loc[pos,'REF'] == max(phase_variant_dict[pos][1], key=phase_variant_dict[pos][1].get):
                    phase_haplo.append('0')
                else:
                    phase_haplo.append('1')
            else:
                phase_haplo.append('0')
            if 2 in phase_variant_dict[pos]:
                if selected_variants.loc[pos,'REF'] == max(phase_variant_dict[pos][2], key=phase_variant_dict[pos][2].get):
                    phase_haplo.append('0')
                else:
                    phase_haplo.append('1')
            else:
                phase_haplo.append('0')
            split_list[0] = '|'.join(phase_haplo)
            selected_variants.loc[pos,'SAMPLE'] = ':'.join(split_list)
                
    
selected_variants.to_csv(vcf_out,mode='a',index=False,sep='\t')

