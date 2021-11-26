import pandas as pd 
import numpy as np
import random
import sys


def sample_data(data, sample_num):
    """
    Sampling data. The number of sampling is sample_num.
    """
    random_index = random.sample(range(0, data.shape[0]), sample_num)
    
    sample_data = data[random_index]

    data = np.delete(data, random_index, 0)

    return data, sample_data


def create_train_test_names(data_file, sample_num):
    """
    Create training and test sets
    """

    data = pd.read_excel(data_file)

    data_0 = []
    data_1 = [] 

    for indexs in data.index:   
        uniprot = data.loc[indexs].values[2]
        position = data.loc[indexs].values[4]

        uniprot_position = uniprot + '_' + str(position)

        functionality = data.loc[indexs].values[5]
        if functionality == 'NP':
            data_0.append(uniprot_position)
        elif functionality in ['FP']:
            data_1.append(uniprot_position) 

    # Randomly sample n samples from Class 0 and class 1 to form a test set
    train_0, test_0 = sample_data(np.array(data_0), sample_num)
    train_1, test_1 = sample_data(np.array(data_1), sample_num)

    # print(train_0.shape, train_1.shape)   
    # print(test_0.shape, test_1.shape)    

    train_data = train_0.tolist() + train_1.tolist()
    test_data = test_0.tolist() + test_1.tolist()

    return train_data, test_data


def generate_1d_data(data_file, train_names, test_names):
    """
    Create train_data_1d and test_data_1d
    """

    data = pd.read_excel(data_file)
    # print(data.columns)

    # Represent the letter of the secondary structure into a number
    ss_mapping = {'-':1,'H':2, 'S':3, 'G':4, 'T':5, 'E':6, 'B':7, 'I':8}
    data['secondary_structure'] = data['secondary_structure'].map(ss_mapping)
    
    data = data.fillna(0)
        
    data = data.replace(np.inf, 1)
    data = data.replace(-np.inf, 0)

    
    target_ptm_type = 'Phosphorylation'
    target_res = ['S', 'T', 'Y']

    train_data, train_label = [], []
    test_data, test_label = [], []

    
    for indexs in data.index:
        res = data.loc[indexs].values[3]
        ptm_types = data.loc[indexs].values[6]
        if res not in target_res:
            continue    
        pdb = data.loc[indexs].values[2]
        position = data.loc[indexs].values[4]
        names = pdb + '_' + str(position)

        if data.loc[indexs].values[5] == 'NP':
            if target_ptm_type in ptm_types:
                content = np.array(data.loc[indexs].values[7:]) 
                if names in train_names:    
                    train_data.append(content)
                    train_label.append(0)
                elif names in test_names:
                    test_data.append(content)
                    test_label.append(0)
        elif data.loc[indexs].values[5] == 'FP':  # if label is 'FP'  
            if target_ptm_type in ptm_types:
                # content = np.hstack((np.array(data1.loc[indexs].values[3]), np.array(data1.loc[indexs].values[7:])))
                content = np.array(data.loc[indexs].values[7:])
                if names in train_names:
                    train_data.append(content)
                    train_label.append(1)
                elif names in test_names:
                    test_data.append(content)
                    test_label.append(1)
    
    train_data, train_label = np.array(train_data), np.array(train_label) 
    test_data, test_label = np.array(test_data), np.array(test_label)   


    np.save('train_data_1d.npy', train_data)
    np.save('train_label_1d.npy', train_label)
    np.save('test_data_1d.npy', test_data)
    np.save('test_label_1d.npy', test_label)      


def get_data(read_file):  
    """
    Change the content of the file into a dictionary
    """

    data = pd.read_excel(read_file)

    ss_mapping = {'-':1,'H':2, 'S':3, 'G':4, 'T':5, 'E':6, 'B':7, 'I':8}
    data['secondary_structure'] = data['secondary_structure'].map(ss_mapping)
    
    data = data.fillna(0)
        
    data = data.replace(np.inf, 1)
    data = data.replace(-np.inf, 0)

    pdb_site_type_feature = {}   # {pdb: {site_type: feature}}
    cur_pdb = ''
    cur_count = 0
    i = 0
    for indexs in data.index:
        pdb = data.loc[indexs].values[2]
        if pdb != cur_pdb:
            cur_pdb = pdb 
            pdb_site_type_feature[cur_pdb] = {}
            cur_count = 0  # first residue
        if pdb == cur_pdb:
            s_type = str(data.loc[indexs].values[3]) + '_' + str(data.loc[indexs].values[4]) + '_' + data.loc[indexs].values[5] + '_' + data.loc[indexs].values[6] 
            cur_count += 1
            site_type = str(cur_count) + '_' + s_type 
            pdb_site_type_feature[cur_pdb][site_type] = np.array(data.loc[indexs].values[7:])
    
    return pdb_site_type_feature


def generate_2d_data(read_dict, window_num, data_names, save_data_file, save_label_file): 
    """
    Read the file, select the ptm site and its left and right window_num sites, 
    and create the training set and test set
    """

    pdb_array = {}    # An array of all the features 
    pdb_site_type = {}  # The corresponding feature of residues
    for pdb in pdb_site_type_feature.keys():
        pdb_array[pdb] = []
        pdb_site_type[pdb] = []
        for site_type in pdb_site_type_feature[pdb].keys():  
            pdb_site_type[pdb].append(site_type)
            pdb_array[pdb].append(pdb_site_type_feature[pdb][site_type])
            
        # break

    data, label = [], []


    res_mapping = {'A':1,'R':2,'N':3,'D':4,'C':5,'G':6,'Q':7,'E':8,'H':9,'I':10,'L':11,
                   'K':12,'M':13,'P':14, 'F':15,'S':16,'T':17,'Y':18,'W':19,'V':20}
    
    target_ptm_type = 'Phosphorylation'
    target_res = ['S', 'T', 'Y']


    types = {'NP':0, 'FP':1}
    window = window_num
    for pdb in pdb_array.keys():
        # print(pdb)
        length = len(pdb_site_type[pdb])
        feature_num = pdb_array[pdb][0].shape[0]  # 65

        for s_t in pdb_site_type[pdb]:  
            site = int(s_t.split('_')[0])  
            res = s_t.split('_')[1]
            position = s_t.split('_')[2]
            sitetype = s_t.split('_')[3]
            ptm_types = s_t.split('_')[4].split(' ')    
            # continue      
            if sitetype in types.keys(): 
                if res not in target_res: 
                    continue
                if (sitetype == 'FP' or sitetype == 'NP') and target_ptm_type not in ptm_types:  # If it is PTM, but not phosphorylation, then skip
                    continue
                pdb_position = pdb + '_' + str(position)
                
                if pdb_position not in data_names:
                    continue
                array = np.zeros((2 * window + 1, feature_num))
                array[window] = pdb_array[pdb][site-1]  # Put the feature vector of the site in the two-dimensional array. List index starts from 0
                cur_window = window  # 6
                cur_site = site # 1
                while cur_site > 1 and cur_window > 0: # Fill the first 6 rows
                    cur_window -= 1
                    cur_site -= 1
                    array[cur_window] = pdb_array[pdb][cur_site-1]
                cur_window = window  # 6
                cur_site = site # 1
                while cur_site < length and cur_window < 2 * window: # Fill the last 6 rows
                    cur_window += 1
                    cur_site += 1
                    array[cur_window] = pdb_array[pdb][cur_site-1]

                data.append(array)
                label.append(types[sitetype])
        

    data, label = np.array(data), np.array(label)   
    
    np.save(save_data_file, data)
    np.save(save_label_file, label)


if __name__ == '__main__':
    data_file = sys.argv[1]
    sample_num = sys.argv[2]
    window_num = sys.argv[3]

    train_names, test_names = create_train_test_names(data_file, sample_num)

    generate_1d_data(data_file, train_names, test_names)

    all_dict = get_data(data_file)  
    generate_2d_data(all_dict, train_names, window_num, 'train_data_2d_all.npy', 'train_label_2d_all.npy')
    generate_2d_data(all_dict, test_names, window_num, 'test_data_2d_all.npy', 'test_label_2d_all.npy')
