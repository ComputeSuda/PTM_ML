import tensorflow as tf
import tensorflow.keras.metrics as metrics
from keras import models
from keras import layers, Input, regularizers 
from keras import optimizers
from keras import activations
from keras.utils import np_utils
from keras import backend as K
from tensorflow.keras import callbacks
import numpy as np 
from sklearn import preprocessing
from sklearn.utils import class_weight
from cal_metric import cal_two_class_metric
import datetime
from keras.models import load_model
import sys
import random


def build_model(dim1, dim2):
    """
    构建模型
    """

    input_x = Input(shape=(dim1, dim2))
    input_xx = layers.BatchNormalization()(input_x)

    x_1 = layers.Conv1D(200, 1, padding='same')(input_xx) # 输出维度200 卷积核1
    x_1 = activations.relu(x_1)

    x_1 = layers.Conv1D(150, 3, padding='same')(x_1)
    x_1 = activations.relu(x_1)

    x_1 = layers.Conv1D(100, 5, padding='same')(x_1)
    x_1 = activations.relu(x_1)

    x_1 = layers.Conv1D(200, 7, padding='same')(x_1)
    x_1 = activations.relu(x_1)
    
    out_x_1 = layers.GlobalAvgPool1D()(x_1)


    x_2 = layers.LSTM(200, return_sequences=True)(input_xx)
    x_2 = layers.LSTM(150, return_sequences=True)(x_2)
    out_x_2 = layers.LSTM(100)(x_2)    

    x = layers.concatenate([out_x_1, out_x_2], axis=-1)
    a = layers.Dense(x.shape[1], activation='softmax')(x) 
    a_probs = a 
    x = layers.Multiply()([x, a_probs])   # (None, 11, 56)

    x = layers.Dense(150, activation='relu')(x)
    x = layers.Dense(100, activation='relu')(x)
    x = layers.Dense(64, activation='relu')(x)
    out_x = layers.Dense(1, activation='sigmoid')(x)

    model = models.Model(input_x, out_x)

    return model


def sample_train_data(data, label):
    """
    对数据进行采样，新的样本中正负样本比例为 1:1
    """
    class_0_index = np.where(label==0)[0]  # 找出0类的索引
    class_1_index = np.where(label==1)[0]

    data_0 = data[class_0_index]  # 类别0的数据
    data_1 = data[class_1_index]

    # 从负样本中采样和正样本一样数量的负样本
    random_class_0 = random.sample(range(0, len(data_0)), len(data_1))

    sample_data = np.vstack([data_0[random_class_0], data_1])
    sample_label = np.zeros(sample_data.shape[0])
    # print(final_data.shape)
    sample_label[:-len(data_1)] = 0
    sample_label[-len(data_1) : ] = 1   

    return sample_data, sample_label


def train_model(train_data_file, train_label_file, models_num):
    """
    在训练集上训练模型
    """

    epochs = 80
    batch_size = 128


    train_data = np.load(train_data_file, allow_pickle=True).astype(float)
    train_label = np.load(train_label_file)

    # 重复models_num次 每次进行 1:1的采样，训练对应的模型
    for model_number in range(models_num):      
        data, label = sample_train_data(train_data, train_label)
        # print(data.shape, label.shape)

        dim1 = data.shape[1]
        dim2 = data.shape[2]
        # sys.exit()

        # print(model_number)
        model = build_model(dim1, dim2)
        model.compile(optimizer=optimizers.RMSprop(lr=0.001), loss='binary_crossentropy', metrics=['acc', metrics.SensitivityAtSpecificity(0.5), metrics.SpecificityAtSensitivity(0.5), metrics.AUC()])
        # print(model.summary())
        # sys.exit()
        
        # callbacks_list = get_callbcak_list(i)
        # cw = class_weight.compute_class_weight('balanced', np.unique(train_label), train_label)
        # cw = dict(enumerate(cw))
        # print(cw)
        # history = model.fit(train_data, train_label, epochs=epochs, batch_size=batch_size, class_weight=cw, verbose=1)
        history = model.fit(data, label, epochs=epochs, batch_size=batch_size, verbose=1)
        # history = history.history
    

        file = str(model_number) + '.h5'
        model.save(file)


def test_model(test_data_file, test_label_file, models_num):
    """
    在测试集上测试模型
    """

    test_data = np.load(test_data_file, allow_pickle=True).astype(float)
    test_label = np.load(test_label_file)

    pred_prob = None
    for model_number in range(models_num):              
        model = load_model(str(model_number) + '.h5', compile=False)

        if model_number == 0:
            pred_prob = model.predict(test_data)
        else:
            pred_prob += model.predict(test_data)

    pred_prob /= models_num
    # print(pred_prob)
    matrix, metric = cal_two_class_metric(test_label, pred_prob)
    print(matrix)
    print(metric)


    with open('result.txt', 'w') as f:
        f.write(str(matrix) + '\n\n')
        for k, v in metric.items():
            f.write(k + ':' + str(v) + '\n')


if __name__ == '__main__':
    train_data_file = sys.argv[1]
    train_label_file = sys.argv[2]
    test_data_file = sys.argv[3]
    test_label_file = sys.argv[4]

    train_model(train_data_file, train_label_file, 5)
    test_model(test_data_file, test_label_file, 5)
