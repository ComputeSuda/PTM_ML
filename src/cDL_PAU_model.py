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
from sklearn.model_selection import KFold
import os 
import sys
import random
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from scipy import interp


def build_model(dim1, dim2):
    """
    Building model
    """

    input_x = Input(shape=(dim1, dim2))
    input_xx = layers.BatchNormalization()(input_x)

    x_1 = layers.Conv1D(200, 1, padding='same')(input_xx) 
    x_1 = activations.relu(x_1)
    x_1 = layers.MaxPooling1D(pool_size=2, strides=1)(x_1)

    x_1 = layers.Conv1D(200, 3, padding='same')(x_1)
    x_1 = activations.relu(x_1)
    x_1 = layers.MaxPooling1D(pool_size=2, strides=1)(x_1)

    x_1 = layers.Conv1D(200, 5, padding='same')(x_1)
    x_1 = activations.relu(x_1)
    x_1 = layers.MaxPooling1D(pool_size=2, strides=1)(x_1)

    x_1 = layers.Conv1D(200, 7, padding='same')(x_1)
    x_1 = activations.relu(x_1)
    x_1 = layers.MaxPooling1D(pool_size=2, strides=1)(x_1)
    
    out_x_1 = layers.GlobalAvgPool1D()(x_1)


    x_2 = layers.LSTM(200, return_sequences=True)(input_xx)
    x_2 = layers.LSTM(100, return_sequences=True)(x_2)
    out_x_2 = layers.LSTM(100)(x_2)    

    x = layers.concatenate([out_x_1, out_x_2], axis=-1)

    x = layers.Dense(256, activation='relu')(x)
    x = layers.Dense(256, activation='relu')(x)
    x = layers.Dense(128, activation='relu')(x)
    x = layers.Dense(64, activation='relu')(x)

    out_x = layers.Dense(1, activation='sigmoid')(x)

    model = models.Model(input_x, out_x)

    return model


def sample_train_data(data, label):
    """
    Sample the data, and the proportion of positive and negative samples in the new sample is 1:1
    """
    class_0_index = np.where(label==0)[0]  # find class 0 indexs
    class_1_index = np.where(label==1)[0]

    data_0 = data[class_0_index]  # class 0 data
    data_1 = data[class_1_index]  # class 1 data

    # Sampling from negative samples and the same number of negative samples
    random_class_0 = random.sample(range(0, len(data_0)), len(data_1))

    sample_data = np.vstack([data_0[random_class_0], data_1])
    sample_label = np.zeros(sample_data.shape[0])
    # print(final_data.shape)
    sample_label[:-len(data_1)] = 0
    sample_label[-len(data_1) : ] = 1   

    return sample_data, sample_label


def combine(cv_scores, metric):
    '''
    The operation of the merger dicts
    '''
    for k,v in metric.items():
        if k not in cv_scores:
            cv_scores[k] = v 
        else:
            cv_scores[k] += v 
    return cv_scores


def get_final_result(cv_scores, k):
    '''
    Get final score
    '''
    for key in cv_scores.keys():
        cv_scores[key] /= k

    return cv_scores


def run_k_fold(train_data, train_label, k):
    '''
    For k-fold cross-validation
    '''
    epochs = 80
    batch_size = 128

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)    

    kfold = KFold(n_splits=k, shuffle=True, random_state=1)
    # Save each folded score
    cv_scores = {}
    cv_count = 0
    dim1 = train_data.shape[1]
    dim2 = train_data.shape[2]

    file_path = './fold_'+ str(k) +'/'
    if not os.path.exists(file_path):
        os.mkdir(file_path)

    for train, test in kfold.split(train_data, train_label):
        cv_count += 1
        print('cv_count:', cv_count)      
        np.save(file_path + 'train_data_'+ str(cv_count) +'.npy', train_data[train])
        np.save(file_path + 'train_label_'+ str(cv_count) +'.npy', train_label[train])
        np.save(file_path + 'test_data_'+ str(cv_count) +'.npy', train_data[test])
        np.save(file_path + 'test_label_'+ str(cv_count) +'.npy', train_label[test])
        # print(train_data[train].shape, train_label[train].shape, np.sum(train_label[train]==0), np.sum(train_label[train]==1))
        # print(train_data[test].shape, train_label[test].shape, np.sum(train_label[test]==0), np.sum(train_label[test]==1))
        
        model = build_model(dim1, dim2)
        model.compile(optimizer=optimizers.RMSprop(lr=0.001), loss='binary_crossentropy', metrics=['acc', metrics.SensitivityAtSpecificity(0.5), metrics.SpecificityAtSensitivity(0.5), metrics.AUC()])

        model.fit(train_data[train], train_label[train], epochs=epochs, batch_size=batch_size, verbose=1)

        pred_prob = model.predict(train_data[test])
        matrix, metric = cal_two_class_metric_capsule(train_label[test], pred_prob)
        # print(metric)

        model.save(file_path + 'fold_' + str(cv_count) + '.h5')
        
        # combine the result
        cv_scores = combine(cv_scores, metric)

        fpr, tpr, thresholds = roc_curve(train_label[test], pred_prob)
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)

    cv_scores = get_final_result(cv_scores, k) 
    with open(file_path + 'fold_' + str(k) + '_result.txt', 'w') as f:
        for k, v in cv_scores.items():
            f.write(k + ':' + str(v) + '\n')     

    # plot roc curve
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)

    plt.plot(mean_fpr, mean_tpr, color=colors[i], label=str(k) + '-fold (AUC = ' + str(mean_auc) + ')', linewidth=2, linestyle=linestyles[i])
    plt.legend(loc = 'lower right', fontsize=14)
    plt.plot([0, 1], [0, 1],'r--', linewidth=2)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.ylabel('True Positive Rate', size=15, fontdict={'weight':'bold'})
    plt.xlabel('False Positive Rate', size=15, fontdict={'weight':'bold'})
    plt.title('ROC', fontsize=15, fontdict={'weight':'bold'})
    plt.xticks(fontsize=15, weight='bold')
    plt.yticks(fontsize=15, weight='bold')
    plt.savefig(file_path + 'fold_' + str(k)+ '_roc.png', dpi=500)
    plt.show()    


def train_model(train_data_file, train_label_file, models_num):
    """
    Training model on training set
    """

    epochs = 80
    batch_size = 128

    train_data = np.load(train_data_file, allow_pickle=True).astype(float)
    train_label = np.load(train_label_file)

    k = 10
    run_k_fold(train_data, train_label, k)

    # The 1:1 sampling was repeated models_num times to train the corresponding model
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
        
        history = model.fit(data, label, epochs=epochs, batch_size=batch_size, verbose=1)
        # history = history.history    

        file = str(model_number) + '.h5'
        model.save(file)


def test_model(test_data_file, test_label_file, models_num):
    """
    Test model on test set
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

    # plor roc
    fpr, tpr, threshold = roc_curve(test_label, pred_prob)
    # roc_auc = str(cut(auc(fpr, tpr),3))
    roc_auc = str(round(auc(fpr, tpr), 3))
    plt.plot(fpr, tpr, '#B38FBC', label = 'FuncPhos(AUC=' + roc_auc + ')', linewidth=2) 
    plt.plot([0, 1], [0, 1],'r--', linewidth=2)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.ylabel('True Positive Rate', size=15, fontdict={'weight':'bold'})
    plt.xlabel('False Positive Rate', size=15, fontdict={'weight':'bold'})
    plt.title('ROC', fontsize=15, fontdict={'weight':'bold'})
    plt.xticks(fontsize=15, weight='bold')
    plt.yticks(fontsize=15, weight='bold')      
    plt.savefig('./FuncPhos.png', dpi=500)
    plt.show() 


if __name__ == '__main__':
    # train_data_file = sys.argv[1]
    # train_label_file = sys.argv[2]
    # test_data_file = sys.argv[3]
    # test_label_file = sys.argv[4]

    # use our datasets    
    train_data_file = '../Datasets/cDL-PAU/Phosphorylation/feature_65/train_data_2d_all.npy'
    train_label_file = '../Datasets/cDL-PAU/Phosphorylation/feature_65/train_label_2d_all.npy'
    test_data_file = '../Datasets/cDL-PAU/Phosphorylation/feature_65/test_data_2d_all.npy'
    test_label_file = '../Datasets/cDL-PAU/Phosphorylation/feature_65/test_data_2d_all.npy'
    

    train_model(train_data_file, train_label_file, 5)
    test_model(test_data_file, test_label_file, 5)
