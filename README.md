This repo is for Functional PTM Predictions
=========================================

- ## 1. Purpose

    To predict Functional-PTM site using deep learning method.

    ## 2.Repo Contents
    - Datasets: All datasets.
    - Evaluation_Indicators: evaluation index functions.
    - fig: plotting ROC curves.
    - src: genarating data, traing or predicting model.

    ## 3. System requirements

    Installation and running has been tested in Ubuntu 18.04.4 LST with python 3.7.7.

    #### Package version

    - python = 3.7.7
    - numpy = 1.18.5
    - Keras = 2.4.3
    - Tensorflow = 2.2.0

    #### You can install the dependent packages by the following commands:

    - pip install python==3.7.7
    - pip install numpy==1.18.5
    - pip install keras==2.4.3
    - pip install tensorflow==2.2.0

    ## 4. Predicting

    - To predict Functional PTM site using all(seq+str+dyn) features, run the following command:

        python ./src/model.py

    ### Training and testing data are provided in the folder of ./Datasets.



This project is licensed under the terms of the GNU General Public License v3.0 license.
