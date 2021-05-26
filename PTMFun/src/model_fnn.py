import tensorflow as tf
import tensorflow.keras.metrics as metrics
from keras import models
from keras import layers, Input, regularizers
from keras import optimizers
from keras import activations
from keras.utils import np_utils
import numpy as np 
from sklearn import preprocessing
from sklearn.utils import class_weight
from cal_metric import cal_two_class_metric
from keras.models import load_model
import sys


test_data = np.load('test_data.npy', allow_pickle=True).astype(float)
test_label = np.load('test_label.npy').astype(float)


def build_model():
	input_x = Input(shape=(1, ))

	x = layers.Dense(16, activation='relu')(input_x)
	x = layers.Dense(8, activation='relu')(x)
	out_x = layers.Dense(1, activation='sigmoid')(x)

	model = models.Model(input_x, out_x)

	return model


def train_model(ite, model_nums):
	epochs = 200
	batch_size = 32

	for i in range(ite):
		train_data = np.load('train_data_' + str(i) + '.npy', allow_pickle=True).astype(float)
		train_label = np.load('train_label_' + str(i) + '.npy').astype(float)

		for j in range(model_nums):
			print(i, j)
			model = build_model()
			model.compile(optimizer=optimizers.RMSprop(lr=0.001), loss='binary_crossentropy', metrics=['acc', metrics.SensitivityAtSpecificity(0.5), metrics.SpecificityAtSensitivity(0.5), metrics.AUC()])

			history = model.fit(train_data, train_label, epochs=epochs, batch_size=batch_size, verbose=0)
			history = history.history

			file = './model/' + str(i) + '/model_fnn_' + str(i) + '_' + str(j) + '.h5'
			model.save(file)


def test_model(ite, model_nums):
	pred_prob = None
	for i in range(ite):
		for j in range(model_nums):
			model = load_model('./model/' + str(i) + '/model_fnn_' + str(i) + '_' + str(j) + '_prl.h5', compile=False)
			pred_prob = model.predict(test_data)			
			matrix, metric = cal_two_class_metric(test_label, pred_prob)
	

if __name__ == '__main__':
	ite = 1

	train_model(ite=6, model_nums=1)
	test_model(ite=6, model_nums=1)



