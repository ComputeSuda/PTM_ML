import numpy as np 
import pandas as pd 
import sys

def generate_1d_data(file):  
	data = pd.read_excel(file)
	
	res_mapping = {'A':1,'R':2,'N':3,'D':4,'C':5,'G':6,'Q':7,'E':8,'H':9,'I':10,'L':11,
				   'K':12,'M':13,'P':14, 'F':15,'S':16,'T':17,'Y':18,'W':19,'V':20}
	data['res'] = data['res'].map(res_mapping)
	
	data = data.fillna(0)
	data = data.replace(np.inf, 1)
	data = data.replace(-np.inf, 0)	
	
	data_0, data_1 = [], []
	for indexs in data.index:		
		if data.loc[indexs].values[5] == 'NP':		
			content = np.array(data.loc[indexs].values[7:])
			data_0.append(content)
		elif data.loc[indexs].values[5] == 'FP':
			content = np.array(data.loc[indexs].values[7:])
			data_1.append(content)

	data_0, data_1 = np.array(data_0), np.array(data_1)  

	return data_0, data_1

	
def sample_data(data, sample_num):
	random_index = random.sample(range(0, data.shape[0]), sample_num)	
	sample_data = data[random_index]
	data = np.delete(data, random_index, 0)
	return data, sample_data


def combine_data(data_0, data_1):
	combine_data = np.vstack([data_0, data_1])
	combine_label = np.zeros(combine_data.shape[0])
	n = data_1.shape[0]

	combine_label[ : -n] = 0
	combine_label[-n : ] = 1

	return combine_data, combine_label

	
def new_sample_1d_data(class_0_data, class_1_data):	
	test_num = int(class_1_data.shape[0] * 0.7)
	class_0_data, sample_0_data = sample_data(class_0_data, test_num)
	class_1_data, sample_1_data = sample_data(class_1_data, test_num)

	test_data, test_label = combine_data(sample_0_data, sample_1_data)
	np.save('test_data.npy', test_data)
	np.save('test_label.npy', test_label)

	nums = class_0_data.shape[0]
	count = 0
	while nums > 0:
		if class_0_data.shape[0]  > class_1_data.shape[0]: 
			sample_num = class_1_data.shape[0]
		else:
			sample_num = class_0_data.shape[0]

		# if (class_0_data.shape[0] - class_1_data.shape[0]) >= class_1_data.shape[0]:
		# 	 sample_num = class_1_data.shape[0]
		# else:
		# 	sample_num = class_0_data.shape[0]

		class_0_data, sample_0_data = sample_data(class_0_data, sample_num)
		train_data, train_label = combine_data(sample_0_data, class_1_data)  
		print(train_data.shape, train_label.shape)
		# print(train_data[1])

		np.save('train_data' + str(count) + '.npy', train_data)
		np.save('train_label' + str(count) + '.npy', train_label)

		count += 1
		nums = class_0_data.shape[0]



if __name__ == '__main__':
	file = sys.argv[1]
	data_0, data_1 = generate_1d_data(file)
	new_sample_1d_data(data_0, data_1)


